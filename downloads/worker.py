"""Bounded HTTP workers with per-resource deadlines and shared origin policy.

Runs in ordinary Python, never imports bpy. The entry point uses os._exit after
closing IPC: live network threads must not keep a disconnected helper alive.
"""
from concurrent.futures import ThreadPoolExecutor
from email.utils import parsedate_to_datetime
import gzip
import http.client
import io
import json
import math
import os
from pathlib import Path
import random
import select
import socket
import ssl
import struct
import sys
import threading
import time
import zlib
from urllib import request, error
from urllib.parse import urlsplit

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from downloads.diagnostics import redact_text, redact_body
from downloads import STALLED_WORKER_EXIT_CODE

WORKER_DEADLINE_GRACE = 5.


class Supervisor:
    """Terminate a helper whose network thread cannot honor its deadline.

    DNS/OS calls cannot always be interrupted by shutting down a socket. This
    thread remains independent of both those calls and blocked IPC writes.
    Finished futures are excluded, even if Blender has not received their bytes.
    """
    def __init__(self, grace):
        self.grace = grace
        self.condition = threading.Condition()
        self.pending = {}
        self.closed = False
        threading.Thread(target=self.run, daemon=True).start()

    def watch(self, future, timeout):
        with self.condition:
            self.pending[future] = time.monotonic() + timeout + self.grace
            self.condition.notify_all()
        future.add_done_callback(self.finished)

    def finished(self, future):
        with self.condition:
            self.pending.pop(future, None)
            self.condition.notify_all()

    def run(self):
        with self.condition:
            while not self.closed:
                # Check done() as well: the callback may be waiting for this lock.
                ends = [end for future, end in self.pending.items() if not future.done()]
                remaining = min(ends) - time.monotonic() if ends else None
                if remaining is not None and remaining <= 0:
                    os._exit(STALLED_WORKER_EXIT_CODE)
                self.condition.wait(remaining)

    def close(self):
        with self.condition:
            self.closed = True
            self.pending.clear()
            self.condition.notify_all()


class Failure(Exception):
    def __init__(self, code, message, **detail):
        self.detail = dict(code=code, message=message, **detail)
        super().__init__(message)


class Deadline:
    def __init__(self, seconds):
        self.end = time.monotonic() + seconds
        self.expired = threading.Event()
        self.lock = threading.Lock()
        self.connections = []
        self.sockets = []
        self.timer = threading.Timer(seconds, self.expire)
        self.timer.daemon = True
        self.timer.start()

    def check(self):
        remaining = self.end - time.monotonic()
        if remaining <= 0 or self.expired.is_set():
            raise Failure('timeout', 'Download deadline exceeded')
        return remaining

    def watch_connection(self, connection):
        with self.lock:
            self.connections.append(connection)
        self.check()

    def watch_socket(self, sock):
        with self.lock:
            self.sockets.append(sock)
        self.check()

    def expire(self):
        self.expired.set()
        # shutdown interrupts reads even if HTTPResponse retains a makefile
        # reference after HTTPConnection has cleared its own socket reference.
        with self.lock:
            sockets = self.sockets + [c.sock for c in self.connections if c.sock is not None]
        for sock in sockets:
            try:
                sock.shutdown(socket.SHUT_RDWR)
            except OSError:
                pass

    def wait(self, seconds):
        self.expired.wait(min(seconds, self.check()))
        self.check()

    def close(self):
        self.timer.cancel()
        for connection in self.connections:
            connection.close()


def retry_delay(exc, attempt):
    value = exc.headers.get('Retry-After') if isinstance(exc, error.HTTPError) else None
    if value:
        try:
            delay = float(value)
        except ValueError:
            try:
                delay = parsedate_to_datetime(value).timestamp() - time.time()
            except (ValueError, TypeError, OverflowError):
                delay = float('nan')
        if math.isfinite(delay):
            return max(0., delay)
    return .25 * 2 ** attempt + random.uniform(0, .1)


class ProviderGate:
    def __init__(self, workers):
        self.workers = workers
        self.condition = threading.Condition()
        self.origins = {}
        self.denied = set()

    def state(self, origin):
        return self.origins.setdefault(origin, dict(active=0, limit=self.workers, until=0., limited=0, status=0))

    def acquire(self, origin, scope, deadline):
        with self.condition:
            state = self.state(origin)
            while True:
                if scope and (origin, scope) in self.denied:
                    raise Failure('authorization_blocked',
                                  'Requests stopped for invalid credentials or an expired session', auth_failed=True)
                remaining = deadline.check()
                delay = state['until'] - time.monotonic()
                if delay >= remaining:
                    raise Failure('provider_cooldown',
                                  'HTTP %s; provider cooldown exceeds this download deadline; retry the import later' % state['status'],
                                  status=state['status'], headers={'Retry-After': str(math.ceil(delay))})
                if delay <= 0 and state['active'] < state['limit']:
                    state['active'] += 1
                    return
                self.condition.wait(min(remaining, max(.01, delay) if delay > 0 else .1))

    def release(self, origin):
        with self.condition:
            self.state(origin)['active'] -= 1
            self.condition.notify_all()

    def cooldown(self, origin, seconds, status):
        with self.condition:
            state = self.state(origin)
            until = time.monotonic() + seconds
            if until > state['until']:
                state['until'], state['status'] = until, status
            if status == 429:
                state['limited'] += 1
                if state['limited'] >= 2:
                    state['limit'] = max(1, state['limit'] // 2)
            self.condition.notify_all()

    def deny(self, origin, scope):
        if scope:
            with self.condition:
                self.denied.add((origin, scope))
                self.condition.notify_all()


def make_opener(deadline, verify_tls):
    class HTTPConnection(http.client.HTTPConnection):
        def connect(self):
            deadline.watch_connection(self)
            self.timeout = min(30., deadline.check())
            super().connect()
            deadline.watch_socket(self.sock)

    class HTTPSConnection(http.client.HTTPSConnection):
        def connect(self):
            deadline.watch_connection(self)
            self.timeout = min(30., deadline.check())
            super().connect()
            deadline.watch_socket(self.sock)

    class HTTPHandler(request.HTTPHandler):
        def http_open(self, req):
            return self.do_open(HTTPConnection, req)

    class HTTPSHandler(request.HTTPSHandler):
        def https_open(self, req):
            return self.do_open(HTTPSConnection, req, context=self._context)

    ctx = ssl.create_default_context() if verify_tls else ssl._create_unverified_context()
    return request.build_opener(HTTPHandler(), HTTPSHandler(context=ctx))


def read_body(response, limit, deadline, *, complete=True):
    content = bytearray()
    # read1 permits deadline checks even if a peer continuously trickles data.
    while len(content) <= limit:
        deadline.check()
        try:
            chunk = response.read1(min(65536, limit + 1 - len(content)))
        except http.client.IncompleteRead:
            raise Failure('incomplete_response', 'Incomplete HTTP response') from None
        deadline.check()
        if not chunk:
            break
        content.extend(chunk)
    if complete:
        if len(content) > limit:
            raise Failure('response_too_large', 'Response exceeds the %s byte size limit' % limit)
    expected = response.headers.get('Content-Length')
    if len(content) <= limit and expected is not None and len(content) < int(expected):
        raise Failure('incomplete_response', 'Incomplete HTTP response (%s of %s bytes)' % (len(content), expected))
    return bytes(content[:limit])


def auth_failure(job, status, body):
    if not job.get('auth_scope'):
        return False
    if status == 401:
        return True
    if status not in (400, 403):
        return False
    try:
        data = json.loads(body)
    except (ValueError, TypeError, RecursionError):
        return False
    def reasons(value):
        if isinstance(value, dict):
            for key, item in value.items():
                if key in ('reason', 'status') and isinstance(item, str):
                    yield item
                yield from reasons(item)
        elif isinstance(value, list):
            for item in value:
                yield from reasons(item)
    try:
        return bool(set(reasons(data)) & set(job.get('auth_error_reasons', ())))
    except RecursionError:
        return False


def perform(job, max_bytes, gate, deadline, attempts):
    url = job['url']
    parsed = urlsplit(url)
    if parsed.scheme not in ('http', 'https') or not parsed.netloc:
        raise Failure('invalid_request', 'Downloads require an HTTP or HTTPS URL')
    origin = (parsed.scheme.lower(), parsed.hostname.lower(), parsed.port or (443 if parsed.scheme == 'https' else 80))
    headers = {'User-Agent': 'Blosm', **job['headers']}
    opener = make_opener(deadline, job['verify_tls'])
    for attempt in range(3):
        gate.acquire(origin, job.get('auth_scope', ''), deadline)
        delay, shared = 0., False
        try:
            attempts[0] = attempt + 1
            try:
                with opener.open(request.Request(url, headers=headers), timeout=min(30., deadline.check())) as response:
                    content = read_body(response, max_bytes, deadline)
                break
            except Failure as exc:
                # Restart the complete GET; never append a new attempt to partial
                # bytes. The same deadline and three-attempt budget still apply.
                if exc.detail['code'] != 'incomplete_response' or attempt == 2:
                    raise
                delay = retry_delay(exc, attempt)
            except error.HTTPError as exc:
                status = exc.code
                delay = retry_delay(exc, attempt)
                shared = status in (429, 503)
                if shared:
                    gate.cooldown(origin, delay, status)
                body, body_error = '', ''
                try:
                    body = read_body(exc, 65536, deadline, complete=False).decode('utf-8', errors='replace')
                except Exception:
                    body_error = '; error response body could not be read completely'
                finally:
                    exc.close()
                denied = auth_failure(job, status, body)
                if denied:
                    gate.deny(origin, job['auth_scope'])
                if status not in (429, 500, 502, 503, 504) or attempt == 2 or denied:
                    raise Failure('http_error', 'HTTP %s%s' % (status, body_error), status=status,
                                  headers={'Retry-After': exc.headers.get('Retry-After', '')},
                                  body=body, auth_failed=denied)
            except (error.URLError, OSError, http.client.HTTPException) as exc:
                deadline.check()
                reason = exc.reason if isinstance(exc, error.URLError) else exc
                if isinstance(reason, ssl.SSLCertVerificationError):
                    raise Failure('tls_error', 'TLS certificate verification failed') from None
                if attempt == 2:
                    code = 'timeout' if isinstance(reason, TimeoutError) else 'connection_error'
                    raise Failure(code, 'Network operation timed out' if code == 'timeout' else 'Connection failed') from None
                delay = retry_delay(exc, attempt)
        finally:
            gate.release(origin)
        if not shared:
            deadline.wait(delay)
    if job['decode_gzip'] and content[:2] == b'\x1f\x8b':
        try:
            with gzip.GzipFile(fileobj=io.BytesIO(content)) as stream:
                content = stream.read(max_bytes + 1)
        except (OSError, EOFError, ValueError, zlib.error):
            raise Failure('invalid_gzip', 'Invalid or incomplete gzip response') from None
        deadline.check()
        if len(content) > max_bytes:
            raise Failure('response_too_large', 'Decoded response exceeds the %s byte size limit' % max_bytes)
    return content


def result(job, max_bytes, gate):
    attempts = [0]
    deadline = Deadline(job['timeout'])
    try:
        return b'B', perform(job, max_bytes, gate, deadline, attempts)
    except Failure as exc:
        detail = exc.detail
    except Exception as exc:
        detail = dict(code='download_error', message='Unexpected download error: ' + type(exc).__name__)
    finally:
        deadline.close()
    if deadline.expired.is_set() and detail.get('status') is None:
        detail.update(code='timeout', message='Download deadline exceeded')
    detail['attempts'] = attempts[0]
    detail['message'] = redact_text(detail['message'], job['url'], job['headers'])
    if 'body' in detail:
        detail['body'] = redact_body(detail['body'], job['url'], job['headers'])
    return b'E', json.dumps(detail).encode()


def read_exact(connection, size):
    chunks = bytearray()
    while len(chunks) < size:
        chunk = connection.recv(size - len(chunks))
        if not chunk:
            raise EOFError
        chunks.extend(chunk)
    return chunks


def main():
    workers, max_bytes = int(sys.argv[3]), int(sys.argv[4])
    pool = ThreadPoolExecutor(max_workers=workers)
    gate = ProviderGate(workers)
    supervisor = Supervisor(WORKER_DEADLINE_GRACE)
    try:
        with socket.create_connection(('127.0.0.1', int(sys.argv[1])), timeout=15) as connection:
            connection.settimeout(None)
            connection.sendall(sys.argv[2].encode('ascii'))
            pending = {}
            heartbeat = time.monotonic()
            while True:
                if select.select([connection], [], [], .005)[0]:
                    size, = struct.unpack('!I', read_exact(connection, 4))
                    if size > 65536 or len(pending) >= workers:
                        raise ValueError('Invalid helper request')
                    job = json.loads(read_exact(connection, size))
                    future = pool.submit(result, job, max_bytes, gate)
                    supervisor.watch(future, job['timeout'])
                    pending[future] = job['id']
                for future in list(pending):
                    if future.done():
                        identity = pending.pop(future)
                        kind, content = future.result()
                        connection.sendall(struct.pack('!cQQ', kind, identity, len(content)))
                        connection.sendall(content)
                        del content, future
                if time.monotonic() >= heartbeat:
                    connection.sendall(struct.pack('!cQQ', b'H', 0, 0))
                    heartbeat = time.monotonic() + 1.
    finally:
        supervisor.close()
        pool.shutdown(wait=False, cancel_futures=True)


def run():
    status = 0
    try:
        main()
    except (EOFError, ConnectionError):
        pass
    except BaseException:
        status = 1
    finally:
        # Python's normal shutdown joins executor threads even with wait=False.
        # IPC is closed; no files need flushing. The OS releases all network I/O.
        os._exit(status)


if __name__ == '__main__':
    run()
