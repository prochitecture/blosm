"""Local scheduling and read-only POST policy contracts; ordinary Python."""
from pathlib import Path
import sys
import unittest
from unittest.mock import Mock, patch
import io
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
sys.path.insert(0, str(Path(__file__).resolve().parent))
from release_target import ADDON_ROOT
sys.path.insert(0, str(ADDON_ROOT))
from progressive.core import Job, LazyTransport, YieldToUI
from downloads import Request
from downloads.worker import perform
from downloads.transport import DownloadTransport
from urllib.error import HTTPError


class OsmSchedulerChecks(unittest.TestCase):
    def test_post_roundtrip_through_download_helper(self):
        class Handler(BaseHTTPRequestHandler):
            def log_message(self, *_):
                pass
            def do_POST(self):
                data = self.rfile.read(int(self.headers['Content-Length']))
                self.send_response(200)
                self.send_header('Content-Length', str(len(data)))
                self.end_headers()
                self.wfile.write(data)
        server = ThreadingHTTPServer(('127.0.0.1', 0), Handler)
        thread = threading.Thread(target=server.serve_forever, daemon=True)
        thread.start()
        transport = DownloadTransport(max_bytes=4096)
        try:
            transport.submit(Request('http://127.0.0.1:%s/query' % server.server_port,
                                     method='POST', body='way(42);out;', retry_safe=True))
            reply = None
            deadline = time.monotonic()+10
            while reply is None and time.monotonic() < deadline:
                reply = transport.poll()
                time.sleep(.005)
            self.assertIsNotNone(reply)
            self.assertEqual(reply[1], b'way(42);out;')
            self.assertIsNone(reply[2])
        finally:
            transport.close()
            server.shutdown()
            server.server_close()
            thread.join(timeout=2)

    def test_cpu_steps_commit_boundaries_and_unknown_commands(self):
        factory = Mock(side_effect=AssertionError('HTTP helper created for local work'))
        visited = []
        def producer():
            visited.append('prepare')
            yield None
            visited.append('commit')
            yield YieldToUI()
            visited.append('finalize')
        job = Job(producer(), LazyTransport(factory))
        job.step()
        self.assertEqual(visited, ['prepare', 'commit'])
        job.step()
        self.assertTrue(job.done)
        self.assertIsNone(job.error)
        self.assertEqual(visited[-1], 'finalize')
        factory.assert_not_called()
        def invalid():
            yield object()
        job = Job(invalid(), LazyTransport(factory))
        job.step()
        self.assertIsInstance(job.error, TypeError)

    def test_post_query_body_and_safe_retry(self):
        from dataclasses import asdict
        for retry_safe, expected_attempts in ((False, 1), (True, 2)):
            attempts = [0]
            gate, deadline = Mock(), Mock()
            deadline.check.return_value = 90
            error = HTTPError('http://example.test/query', 503, 'busy', {}, io.BytesIO(b''))
            response = Mock()
            response.__enter__ = Mock(return_value=response)
            response.__exit__ = Mock(return_value=False)
            opener = Mock()
            opener.open.side_effect = [error, response]
            job = asdict(Request('http://example.test/query', method='POST', body='way(1);out;', retry_safe=retry_safe))
            with patch('downloads.worker.make_opener', return_value=opener), \
                    patch('downloads.worker.read_body', side_effect=lambda r, *a, **kw: b'ok' if r is response else b''):
                if retry_safe:
                    self.assertEqual(perform(job, 100, gate, deadline, attempts), b'ok')
                else:
                    with self.assertRaises(Exception):
                        perform(job, 100, gate, deadline, attempts)
            self.assertEqual(opener.open.call_count, expected_attempts)
            request = opener.open.call_args.args[0]
            self.assertEqual(request.get_method(), 'POST')
            self.assertEqual(request.data, b'way(1);out;')


if __name__ == '__main__':
    unittest.main()
