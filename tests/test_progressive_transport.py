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
from collections import Counter
from urllib.error import HTTPError

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from progressive.core import Job, Download, Drain
from downloads import Request
from downloads.transport import DownloadTransport, reap

PAYLOAD = bytes(range(256)) * 16384


class Handler(BaseHTTPRequestHandler):
    def log_message(self, *_):
        pass

    def do_GET(self):
        with self.server.guard:
            self.server.active += 1
            self.server.peak = max(self.server.peak, self.server.active)
            self.server.requests[self.path] += 1
        try:
            self.serve_payload()
        finally:
            with self.server.guard:
                self.server.active -= 1

    def serve_payload(self):
        if self.path.startswith('/work/'):
            index = int(self.path.rsplit('/', 1)[1])
            time.sleep(.25 if index == 0 else .04)
            payload = str(index).encode() * 1024
        else:
            payload = gzip.compress(PAYLOAD) if self.path == '/gzip' else PAYLOAD
        if self.path == '/retry' and self.server.requests[self.path] == 1:
            self.send_response(429)
            self.send_header('Retry-After', '0')
            self.end_headers()
            return
        if self.path == '/missing':
            self.send_error(404)
            return
        if self.path == '/slow':
            time.sleep(2)
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
        cls.server.guard = threading.Lock()
        cls.server.active = cls.server.peak = 0
        cls.server.requests = Counter()
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
                identity = transport.submit(Request(self.url + path, decode_gzip=True))
                deadline = time.monotonic() + 10
                while time.monotonic() < deadline:
                    reply = transport.poll()
                    if reply is not None:
                        break
                    time.sleep(.005)
                else:
                    self.fail('No helper reply')
                received, data, error = reply
                self.assertEqual(identity, received)
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
        for _ in range(4):
            transport.submit(self.url + '/slow')
        for _ in range(20):
            self.assertIsNone(transport.poll())
            time.sleep(.01)
        started = time.monotonic()
        transport.close()
        self.assertLess(time.monotonic() - started, .3)
        process.wait(timeout=5)
        reap()
        # The server threads finish after client cancellation; do not let them
        # inflate the following test's measured peak request count.
        until = time.monotonic() + 3
        while self.server.active and time.monotonic() < until:
            time.sleep(.01)
        transport = DownloadTransport()
        try:
            transport.process.kill()
            transport.process.wait(timeout=5)
            with self.assertRaisesRegex(RuntimeError, 'stopped unexpectedly'):
                transport.poll()
        finally:
            transport.close()

    def test_cancel_discards_buffered_results(self):
        transport = DownloadTransport(max_requests=2)
        process = transport.process
        consumed = []
        def steps():
            for index in range(6):
                yield Download(self.url + '/work/%s' % index,
                               lambda data, error: consumed.append(data))
            yield Drain()
        job = Job(steps(), transport)
        try:
            job.step(budget=.1)
            until = time.monotonic() + 5
            while len(job.replies) < 2 and time.monotonic() < until:
                job.poll()
                time.sleep(.005)
            self.assertEqual(len(job.replies), 2)
            job.cancel()
            job.step()
            self.assertEqual(consumed, [])
            self.assertFalse(job.pending or job.replies)
            self.assertTrue(job.done and job.cancelled)
        finally:
            job.cancel()
            process.wait(timeout=5)

    def test_job_errors_cancel_and_cooperative_bound(self):
        class Transport:
            max_requests = 1
            closed = False
            def submit(self, request):
                self.request = request
                return 0
            def poll(self):
                return 0, None, ValueError('network')
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
        self.assertEqual(events[0], 'handled')
        self.assertLessEqual(events.count('step'), 8)
        self.assertFalse(job.done)
        job.cancel()
        self.assertTrue(job.done and job.cancelled and transport.closed)
        self.assertEqual(events[-1], 'closed')

    def test_concurrency_out_of_order_and_pause_backpressure(self):
        transport = DownloadTransport(max_requests=4, max_bytes=8192, max_buffer_bytes=3 * 8192)
        self.assertEqual(transport.max_requests, 3)
        results, events = [], []
        with self.server.guard:
            self.server.peak = 0
        def consume(index, data, error):
            self.assertIsNone(error)
            self.assertEqual(data, str(index).encode() * 1024)
            results.append(index)
        def steps():
            for index in range(9):
                yield Download(self.url + '/work/%s' % index,
                               lambda data, error, index=index: consume(index, data, error))
            yield Drain()
            events.append('finalized')
            return 9
        job = Job(steps(), transport)
        process = transport.process
        try:
            job.step(budget=.1)
            # Polling while paused must not consume or submit more jobs.
            until = time.monotonic() + 5
            while len(job.replies) != 3 and time.monotonic() < until:
                job.poll()
                time.sleep(.005)
            self.assertEqual(len(job.replies), 3)
            self.assertEqual(results, [])
            self.assertEqual(transport.next_id, 3)
            while not job.done and time.monotonic() < until:
                job.step()
                time.sleep(.001)
            self.assertTrue(job.done)
            self.assertIsNone(job.error)
            self.assertEqual(sorted(results), list(range(9)))
            self.assertNotEqual(results[0], 0)
            self.assertEqual(events, ['finalized'])
            self.assertEqual(job.result, 9)
            self.assertGreater(self.server.peak, 1)
            self.assertLessEqual(self.server.peak, 3)
        finally:
            job.cancel()
            process.wait(timeout=5)

    def test_one_worker_retries_raw_gzip_and_size_limits(self):
        for path, decode, limit in (('/work/1', False, 8192), ('/retry', False, len(PAYLOAD)),
                                    ('/gzip', False, len(PAYLOAD)), ('/gzip', True, 8192)):
            transport = DownloadTransport(max_requests=1, max_bytes=limit)
            received = []
            def steps():
                try:
                    received.append((yield Request(self.url + path, decode_gzip=decode)))
                except RuntimeError as exc:
                    received.append(exc)
            job = Job(steps(), transport)
            try:
                until = time.monotonic() + 10
                while not job.done and time.monotonic() < until:
                    job.step()
                    time.sleep(.002)
                self.assertTrue(job.done)
                self.assertIsNone(job.error)
                if path == '/gzip' and decode:
                    self.assertIsInstance(received[0], RuntimeError)
                elif path == '/gzip':
                    self.assertEqual(gzip.decompress(received[0]), PAYLOAD)
                elif path == '/retry':
                    self.assertEqual(self.server.requests[path], 2)
                    self.assertEqual(received[0], PAYLOAD)
            finally:
                job.cancel()

    def test_refill_overlaps_slow_consumer(self):
        transport = DownloadTransport(max_requests=2)
        seen = []
        with self.server.guard:
            self.server.requests.clear()
        def consume(index, data, error):
            self.assertIsNone(error)
            time.sleep(.15)  # Stand-in for a synchronous Blender import.
            seen.append(index)
            if len(seen) == 2:
                self.assertIn('/work/13', self.server.requests,
                              'The next download must start before ready results are all imported')
        def steps():
            for index in range(11, 17):
                yield Download(self.url + '/work/%s' % index,
                               lambda data, error, index=index: consume(index, data, error))
            yield Drain()
        job = Job(steps(), transport)
        try:
            until = time.monotonic() + 10
            while not job.done and time.monotonic() < until:
                job.step()
                time.sleep(.001)
            self.assertTrue(job.done)
            self.assertIsNone(job.error, job.error)
            self.assertEqual(sorted(seen), list(range(11, 17)))
        finally:
            job.cancel()


if __name__ == '__main__':
    unittest.main()
