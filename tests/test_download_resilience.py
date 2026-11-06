"""Standalone fault-injection tests for shared download resilience."""
from collections import Counter
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json
from pathlib import Path
import sys
import threading
import time
import unittest
from unittest.mock import patch
from urllib.error import HTTPError

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from downloads import Request, DownloadError
from downloads.transport import DownloadTransport
from downloads import transport as transport_module
from downloads.diagnostics import redact_text, redact_url, redact_body


class Handler(BaseHTTPRequestHandler):
    def log_message(self, *_):
        pass

    def do_GET(self):
        path = self.path.split('?', 1)[0]
        with self.server.guard:
            self.server.requests[path] += 1
            number = self.server.requests[path]
            self.server.events.append((path, time.monotonic()))
            self.server.active += 1
            self.server.peak = max(self.server.peak, self.server.active)
        try:
            self.respond(path, number)
        except (ConnectionError, OSError):
            pass
        finally:
            with self.server.guard:
                self.server.active -= 1

    def respond(self, path, number):
        if path == '/hang_headers':
            time.sleep(3)
        if path == '/trickle_headers':
            self.wfile.write(b'HTTP/1.1 200 OK\r\nX-Trickle: ')
            for _ in range(100):
                self.wfile.write(b'a')
                time.sleep(.03)
            return
        if path == '/quota' and number == 1:
            self.send_response(429)
            self.send_header('Retry-After', '1')
            self.send_header('Content-Length', '0')
            self.end_headers()
            return
        if path == '/cooldown_long':
            self.send_response(429)
            self.send_header('Retry-After', '120')
            self.send_header('Content-Length', '0')
            self.end_headers()
            return
        if path == '/quota_twice' and number <= 2:
            self.send_response(429)
            self.send_header('Retry-After', '0')
            self.send_header('Content-Length', '0')
            self.end_headers()
            return
        if path == '/busy':
            time.sleep(.12)
        status = 200
        payload = b'OK'
        if path in ('/error_body_timeout', '/denied', '/isolated', '/confirmed', '/expired', '/bad_request', '/broken_error'):
            status = 400 if path in ('/expired', '/bad_request') else (401 if path == '/denied' else 403)
            payload = json.dumps({'error': {'reason': 'INVALID_KEY' if path in ('/confirmed', '/expired') else 'PERMISSION_DENIED',
                                  'message': 'Failed URL https://host/tiles?key=topsecret&session=private-session; Bearer authsecret'}}).encode()
        if path == '/deep_error':
            status, payload = 403, b'[' * 2000 + b'0' + b']' * 2000
        if path == '/repeated_secret':
            status, payload = 403, b'x' * 65000
        if path == '/large':
            payload = b'x' * 1024
        if path == '/bad_gzip':
            payload = b'\x1f\x8b' + b'not a gzip stream'
        if path == '/trickle_body':
            payload = b'x' * 100
        if path == '/large_ok':
            payload = b'x' * (4 * 1024 * 1024)
        if path == '/retry_chunked':
            self.send_response(status)
            self.send_header('Transfer-Encoding', 'chunked')
            self.end_headers()
            self.wfile.write(b'2\r\nO' if number == 1 else b'2\r\nOK\r\n0\r\n\r\n')
            return
        self.send_response(status)
        self.send_header('Content-Length', str(len(payload)))
        self.end_headers()
        if path in ('/broken_error', '/truncated') or (path == '/retry_truncated' and number == 1):
            self.wfile.write(payload[:1])
            return
        if path == '/error_body_timeout':
            self.wfile.write(payload[:1])
            time.sleep(3)
            return
        if path == '/trickle_body':
            for value in payload:
                self.wfile.write(bytes([value]))
                time.sleep(.03)
            return
        self.wfile.write(payload)


class ResilienceChecks(unittest.TestCase):
    def setUp(self):
        self.servers = []
        self.transports = []
        for _ in range(2):
            server = ThreadingHTTPServer(('127.0.0.1', 0), Handler)
            server.guard = threading.Lock()
            server.requests = Counter()
            server.events = []
            server.active = server.peak = 0
            thread = threading.Thread(target=server.serve_forever, daemon=True)
            thread.start()
            self.servers.append((server, thread))
        self.url, self.other_url = ['http://127.0.0.1:%s' % s.server_port for s, _ in self.servers]

    def tearDown(self):
        for transport in self.transports:
            process = transport.process
            transport.close()
            if process:
                process.wait(timeout=3)
        for server, thread in self.servers:
            server.shutdown()
            server.server_close()
            thread.join()

    def transport(self, **kwargs):
        transport = DownloadTransport(**kwargs)
        self.transports.append(transport)
        return transport

    def supervised_transport(self, mode):
        launch = transport_module.subprocess.Popen
        def fixture(args, **kwargs):
            args = list(args)
            args[3:4] = [str(Path(__file__).with_name('stalled_download_worker.py')), mode]
            return launch(args, **kwargs)
        with patch.object(transport_module.subprocess, 'Popen', side_effect=fixture):
            return self.transport()

    def replies(self, transport, count=1):
        replies = []
        until = time.monotonic() + 6
        while len(replies) < count and time.monotonic() < until:
            reply = transport.poll()
            if reply is not None:
                replies.append(reply)
            time.sleep(.002)
        self.assertEqual(len(replies), count)
        return replies

    def test_deadlines_isolate_hung_headers_and_slow_streams(self):
        for path in ('/hang_headers', '/trickle_headers', '/trickle_body'):
            with self.subTest(path=path):
                transport = self.transport(max_requests=2)
                slow = transport.submit(Request(self.url + path, timeout=.3))
                good = transport.submit(self.url + '/ok')
                start = time.monotonic()
                replies = {identity: (data, failure) for identity, data, failure in self.replies(transport, 2)}
                self.assertLess(time.monotonic() - start, 1.5)
                self.assertEqual(replies[good], (b'OK', None))
                self.assertEqual(replies[slow][1].code, 'timeout')
                self.assertEqual(replies[slow][1].attempts, 1)
                transport.submit(self.url + '/ok')
                self.assertEqual(self.replies(transport)[0][1:], (b'OK', None))

    def test_blender_delay_does_not_expire_completed_reply(self):
        transport = self.transport()
        while not transport.authenticated:
            transport.poll()
            time.sleep(.002)
        identity = transport.submit(Request(self.url + '/ok', timeout=.3))
        transport.poll()  # Send, then simulate Blender being busy longer than the deadline.
        time.sleep(.6)
        transport.last_poll -= 1000
        transport.silent_poll_seconds = 14.9
        self.assertEqual(self.replies(transport)[0], (identity, b'OK', None))

    def test_shared_cooldown_does_not_delay_other_origin(self):
        transport = self.transport(max_requests=4)
        quota = transport.submit(self.url + '/quota')
        server = self.servers[0][0]
        until = time.monotonic() + 3
        while not server.requests['/quota'] and time.monotonic() < until:
            transport.poll()
            time.sleep(.002)
        time.sleep(.1)  # Let the worker observe the first 429 before submitting siblings.
        same = transport.submit(self.url + '/ok')
        other = transport.submit(self.other_url + '/ok')
        replies = self.replies(transport, 3)
        self.assertEqual({r[0] for r in replies}, {quota, same, other})
        self.assertTrue(all(r[1:] == (b'OK', None) for r in replies))
        events = server.events
        first = next(t for path, t in events if path == '/quota')
        later = next(t for path, t in events if path == '/ok')
        other_time = self.servers[1][0].events[0][1]
        self.assertGreaterEqual(later - first, .9)
        self.assertLess(other_time - first, .8)
        self.assertEqual(server.requests['/quota'], 2)

    def test_repeated_rate_limits_reduce_concurrency(self):
        transport = self.transport(max_requests=4)
        transport.submit(self.url + '/quota_twice')
        self.assertEqual(self.replies(transport)[0][1:], (b'OK', None))
        server = self.servers[0][0]
        with server.guard:
            server.peak = 0
        for _ in range(4):
            transport.submit(self.url + '/busy')
        self.replies(transport, 4)
        self.assertLessEqual(server.peak, 2)
        self.assertGreater(server.peak, 1)

    def test_long_retry_after_is_not_retried_early(self):
        transport = self.transport()
        transport.submit(Request(self.url + '/cooldown_long', timeout=.3))
        failure = self.replies(transport)[0][2]
        self.assertIsInstance(failure, HTTPError)
        self.assertEqual(failure.code, 429)
        self.assertEqual(failure.download_code, 'provider_cooldown')
        self.assertEqual(failure.attempts, 1)
        transport.submit(Request(self.url + '/ok', timeout=.3))
        failure = self.replies(transport)[0][2]
        self.assertEqual(failure.attempts, 0)
        self.assertEqual(self.servers[0][0].requests['/cooldown_long'], 1)
        self.assertEqual(self.servers[0][0].requests['/ok'], 0)
        transport.submit(self.other_url + '/ok')
        self.assertEqual(self.replies(transport)[0][1:], (b'OK', None))

    def test_confirmed_auth_failure_blocks_only_its_scope(self):
        transport = self.transport()
        for path, scope, denied in (('/denied', 'a', True), ('/isolated', 'b', False),
                                    ('/confirmed', 'c', True), ('/expired', 'd', True), ('/bad_request', 'e', False)):
            transport.submit(Request(self.url + path + '?key=topsecret&session=private-session',
                                     headers={'Authorization': 'Bearer authsecret'}, auth_scope=scope,
                                     auth_error_reasons=('INVALID_KEY',)))
            failure = self.replies(transport)[0][2]
            self.assertIsInstance(failure, HTTPError)
            self.assertEqual(failure.attempts, 1)
            self.assertEqual(failure.auth_failed, denied)
            text = str(failure) + failure.url + failure.read().decode()
            for secret in ('topsecret', 'private-session', 'authsecret'):
                self.assertNotIn(secret, text)
            before = self.servers[0][0].requests['/ok']
            transport.submit(Request(self.url + '/ok', auth_scope=scope))
            _, data, failure = self.replies(transport)[0]
            if denied:
                self.assertEqual(failure.code, 'authorization_blocked')
                self.assertEqual(failure.attempts, 0)
                self.assertEqual(self.servers[0][0].requests['/ok'], before)
            else:
                self.assertEqual((data, failure), (b'OK', None))
        transport.submit(Request(self.url + '/ok', auth_scope='different-key'))
        self.assertEqual(self.replies(transport)[0][1:], (b'OK', None))

    def test_unusual_error_bodies_stay_bounded_and_isolated(self):
        transport = self.transport()
        for path in ('/deep_error', '/repeated_secret?key=x'):
            transport.submit(self.url + path)
            failure = self.replies(transport)[0][2]
            self.assertIsInstance(failure, HTTPError)
            self.assertEqual(failure.code, 403)
            self.assertLessEqual(len(failure.read()), 65536)
            transport.submit(self.url + '/ok')
            self.assertEqual(self.replies(transport)[0][1:], (b'OK', None))

    def test_error_body_faults_do_not_kill_helper(self):
        transport = self.transport()
        for path in ('/broken_error', '/error_body_timeout'):
            transport.submit(Request(self.url + path, timeout=.3))
            failure = self.replies(transport)[0][2]
            self.assertIsInstance(failure, HTTPError)
            self.assertEqual(failure.code, 403)
            self.assertIn('body could not be read completely', str(failure))
            transport.submit(self.url + '/ok')
            self.assertEqual(self.replies(transport)[0][1:], (b'OK', None))

    def test_specific_diagnostics(self):
        transport = self.transport(max_bytes=64)
        for path, code in (('/large', 'response_too_large'), ('/bad_gzip', 'invalid_gzip'),
                           ('/truncated', 'incomplete_response')):
            transport.submit(Request(self.url + path, decode_gzip=True))
            failure = self.replies(transport)[0][2]
            self.assertIsInstance(failure, DownloadError)
            self.assertEqual(failure.code, code)
            self.assertEqual(failure.attempts, 3 if code == 'incomplete_response' else 1)

    def test_incomplete_transfer_retries_from_scratch(self):
        transport = self.transport()
        for path in ('/retry_truncated', '/retry_chunked'):
            transport.submit(self.url + path)
            self.assertEqual(self.replies(transport)[0][1:], (b'OK', None))
            self.assertEqual(self.servers[0][0].requests[path], 2)

    def test_incomplete_retries_share_the_original_deadline(self):
        transport = self.transport()
        start = time.monotonic()
        transport.submit(Request(self.url + '/truncated', timeout=.1))
        failure = self.replies(transport)[0][2]
        self.assertEqual(failure.code, 'timeout')
        self.assertEqual(failure.attempts, 1)
        self.assertLess(time.monotonic() - start, 1.5)
        transport.submit(self.url + '/ok')
        self.assertEqual(self.replies(transport)[0][1:], (b'OK', None))

    def test_supervisor_stops_stalled_dns_even_with_heartbeats(self):
        transport = self.supervised_transport('dns')
        process = transport.process
        transport.submit(self.url + '/ok')
        self.assertEqual(self.replies(transport)[0][1:], (b'OK', None))
        transport.submit(Request('http://stalled.invalid/tile.glb?key=private-key', timeout=1.1))
        start = time.monotonic()
        with self.assertRaisesRegex(TimeoutError, 'worker remained blocked past its deadline') as caught:
            self.replies(transport)
        self.assertLess(time.monotonic() - start, 3.)
        self.assertNotIn('private-key', str(caught.exception))
        process.wait(timeout=1)
        self.assertEqual(process.returncode, transport_module.STALLED_WORKER_EXIT_CODE)

    def test_supervisor_ignores_completed_bytes_waiting_for_blender(self):
        transport = self.supervised_transport('normal')
        # Establish IPC first; a 4 MiB result then fills the local socket while
        # the parent is paused beyond both the resource deadline and grace.
        transport.submit(self.url + '/ok')
        self.replies(transport)
        transport.submit(Request(self.url + '/large_ok', timeout=.3))
        transport.poll()
        time.sleep(1.)
        _, content, failure = self.replies(transport)[0]
        self.assertIsNone(failure)
        self.assertEqual(content, b'x' * (4 * 1024 * 1024))
        transport.submit(self.url + '/ok')
        self.assertEqual(self.replies(transport)[0][1:], (b'OK', None))

    def test_windows_reset_preserves_stalled_worker_diagnosis(self):
        transport = self.transport()
        transport.submit(self.url + '/ok')
        self.replies(transport)
        # Simulate the ordering observed on Windows: socket reset is visible one
        # poll before the process termination code becomes available.
        with patch.object(transport, 'connection', wraps=transport.connection) as connection, \
                patch.object(transport.process, 'poll', side_effect=[None, transport_module.STALLED_WORKER_EXIT_CODE]):
            connection.recv.side_effect = ConnectionResetError(10054, 'Connection reset')
            self.assertIsNone(transport.poll())
            with self.assertRaisesRegex(TimeoutError, 'worker remained blocked past its deadline'):
                transport.poll()

    def test_helper_exits_on_parent_disconnect_with_active_request(self):
        transport = self.transport()
        process = transport.process
        transport.submit(Request(self.url + '/hang_headers', timeout=30))
        until = time.monotonic() + 3
        while not self.servers[0][0].requests['/hang_headers'] and time.monotonic() < until:
            transport.poll()
            time.sleep(.002)
        self.assertTrue(self.servers[0][0].requests['/hang_headers'])
        transport.connection.close()  # No kill(): the helper must notice EOF itself.
        transport.connection = None
        process.wait(timeout=1.5)
        self.assertEqual(process.returncode, 0)

    def test_redaction_keeps_json_valid(self):
        url = 'https://username:password@host/tile.glb?key=quoted%22secret&session=private-session#fragment'
        text = redact_text('Problem at ' + url, url)
        self.assertIn('host/tile.glb', text)
        self.assertNotIn('username', text)
        self.assertNotIn('password', text)
        self.assertNotIn('private-session', text)
        body = json.dumps({'error': {'message': 'Rejected quoted"secret and Bearer authsecret',
                                    'reason': 'INVALID_KEY', 'access_token': 'unknown-token'}})
        clean = json.loads(redact_body(body, url, {'Authorization': 'Bearer authsecret'}))
        self.assertEqual(clean['error']['reason'], 'INVALID_KEY')
        self.assertNotIn('authsecret', str(clean))
        self.assertNotIn('quoted"secret', str(clean))
        self.assertNotIn('unknown-token', str(clean))


if __name__ == '__main__':
    unittest.main()
