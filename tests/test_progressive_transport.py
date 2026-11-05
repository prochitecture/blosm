"""Run with standalone Python: helper IPC, HTTP and scheduler lifecycle tests.

The HTTP fixture uses threads in this standalone test process, never in Blender.
"""
import gzip
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
import sys
import threading
import time
import unittest
from urllib.error import HTTPError

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from progressive.core import Job
from progressive.transport import DownloadTransport, reap

PAYLOAD = bytes(range(256)) * 16384


class Handler(BaseHTTPRequestHandler):
    def log_message(self, *_):
        pass

    def do_GET(self):
        if self.path == '/missing':
            self.send_error(404)
            return
        if self.path == '/slow':
            time.sleep(2)
        payload = gzip.compress(PAYLOAD) if self.path == '/gzip' else PAYLOAD
        self.send_response(200)
        self.send_header('Content-Length', str(len(payload)))
        self.end_headers()
        try:
            self.wfile.write(payload[:10] if self.path == '/truncated' else payload)
        except (ConnectionError, OSError):
            pass


class TransportChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.server = ThreadingHTTPServer(('127.0.0.1', 0), Handler)
        cls.thread = threading.Thread(target=cls.server.serve_forever, daemon=True)
        cls.thread.start()
        cls.url = 'http://127.0.0.1:%s' % cls.server.server_port

    @classmethod
    def tearDownClass(cls):
        cls.server.shutdown()
        cls.server.server_close()
        cls.thread.join()

    def test_memory_http_gzip_errors_and_reuse(self):
        transport = DownloadTransport()
        try:
            for path in ('/data', '/gzip', '/missing', '/truncated', '/data'):
                transport.submit(self.url + path)
                deadline = time.monotonic() + 10
                while time.monotonic() < deadline:
                    reply = transport.poll()
                    if reply is not None:
                        break
                    time.sleep(.005)
                else:
                    self.fail('No helper reply')
                data, error = reply
                if path == '/missing':
                    self.assertIsInstance(error, HTTPError)
                    self.assertEqual(error.code, 404)
                elif path == '/truncated':
                    self.assertIsInstance(error, RuntimeError)
                else:
                    self.assertIsNone(error)
                    self.assertEqual(data, PAYLOAD)
        finally:
            transport.close()

    def test_cancel_during_download_and_helper_exit(self):
        transport = DownloadTransport()
        process = transport.process
        transport.submit(self.url + '/slow')
        for _ in range(20):
            self.assertIsNone(transport.poll())
            time.sleep(.01)
        started = time.monotonic()
        transport.close()
        self.assertLess(time.monotonic() - started, .3)
        process.wait(timeout=5)
        reap()
        transport = DownloadTransport()
        try:
            transport.process.kill()
            transport.process.wait(timeout=5)
            with self.assertRaisesRegex(RuntimeError, 'stopped unexpectedly'):
                transport.poll()
        finally:
            transport.close()

    def test_job_errors_cancel_and_cooperative_bound(self):
        class Transport:
            closed = False
            def submit(self, request):
                self.request = request
            def poll(self):
                return None, ValueError('network')
            def close(self):
                self.closed = True
        events = []
        def steps():
            try:
                try:
                    yield 'request'
                except ValueError:
                    events.append('handled')
                for _ in range(1000):
                    events.append('step')
                    yield None
                return 42
            finally:
                events.append('closed')
        transport = Transport()
        job = Job(steps(), transport)
        job.step()
        job.step(max_steps=4)
        self.assertEqual(events, ['handled'] + ['step'] * 4)
        self.assertFalse(job.done)
        job.cancel()
        self.assertTrue(job.done and job.cancelled and transport.closed)
        self.assertEqual(events[-1], 'closed')


if __name__ == '__main__':
    unittest.main()
