"""Bounded, memory-only HTTP transport using one isolated worker process.

submit returns an identity; poll returns (identity, bytes, exception) in completion
order. The caller owns interpretation and must bound results retained after poll.
No threads or blocking network reads run in the Blender/caller process.
"""
import atexit
from dataclasses import asdict
import io
import json
import math
from pathlib import Path
import secrets
import socket
import struct
import subprocess
import sys
from time import monotonic
from urllib.error import HTTPError
from . import Request, DownloadError, STALLED_WORKER_EXIT_CODE
from .diagnostics import redact_url
from .access import require_online_access

MAX_BYTES = 256 * 1024 * 1024
_retired = []
_live = set()


def reap():
    _retired[:] = [process for process in _retired if process.poll() is None]
    return bool(_retired)


@atexit.register
def _reap_at_exit():
    for transport in list(_live):
        transport.close()
    for process in _retired:
        process.wait(timeout=5)


def python_executable():
    root = Path(sys.prefix)
    for candidate in (root / 'bin' / 'python.exe', root / 'bin' / 'python3',
                      root / 'bin' / ('python%d.%d' % sys.version_info[:2])):
        if candidate.is_file():
            return str(candidate)
    if Path(sys.executable).name.lower().startswith('python'):
        return sys.executable
    raise RuntimeError("Could not locate Blender's bundled Python for downloading")


class DownloadTransport:
    def __init__(self, max_requests=4, max_bytes=MAX_BYTES, max_buffer_bytes=None):
        require_online_access()
        if not 1 <= max_requests <= 8 or not 1 <= max_bytes <= MAX_BYTES:
            raise ValueError('Invalid download limits')
        # Reserve a worst-case payload for every outstanding request. Thus
        # unknown Content-Length and compressed replies cannot overrun the budget.
        self.max_buffer_bytes = max_buffer_bytes if max_buffer_bytes is not None else max_requests * max_bytes
        self.max_requests = min(max_requests, self.max_buffer_bytes // max_bytes)
        if not self.max_requests:
            raise ValueError('Download memory budget must fit at least one response')
        self.max_bytes = max_bytes
        self.listener = socket.socket()
        self.connection = None
        self.process = None
        self.token = secrets.token_hex(16).encode('ascii')
        self.authenticated = False
        self.incoming = bytearray()
        self.outgoing = bytearray()
        self.pending = {}
        self.next_id = 0
        self.last_poll = monotonic()
        self.silent_poll_seconds = 0.
        self.disconnected = None
        try:
            self.listener.bind(('127.0.0.1', 0))
            self.listener.listen(1)
            self.listener.setblocking(False)
            self.process = subprocess.Popen(
                [python_executable(), '-I', '-B', str(Path(__file__).with_name('worker.py')),
                 str(self.listener.getsockname()[1]), self.token.decode('ascii'),
                 str(self.max_requests), str(max_bytes)],
                stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                creationflags=getattr(subprocess, 'CREATE_NO_WINDOW', 0))
            _live.add(self)
        except Exception:
            self.close()
            raise

    def submit(self, request):
        self._check_online_access()
        if self.process is None:
            raise RuntimeError('Download transport is closed')
        if len(self.pending) >= self.max_requests:
            raise RuntimeError('Download capacity exhausted')
        if isinstance(request, str):
            request = Request(request)
        if request.method not in ('GET', 'POST') or not isinstance(request.body, str):
            raise ValueError('Downloads support GET or POST with a text body')
        if not math.isfinite(request.timeout) or not 0 < request.timeout <= 3600:
            raise ValueError('Download timeout must be between 0 and 3600 seconds')
        identity = self.next_id
        data = json.dumps(dict(asdict(request), id=identity)).encode('utf-8')
        if len(data) > 65536:
            raise ValueError('Download request is too long')
        self.next_id += 1
        self.pending[identity] = request.url
        self.outgoing.extend(struct.pack('!I', len(data)) + data)
        return identity

    def poll(self):
        self._check_online_access()
        if self.process is None:
            raise RuntimeError('Download transport is closed')
        returncode = self.process.poll()
        if returncode == STALLED_WORKER_EXIT_CODE:
            raise TimeoutError('A download worker remained blocked past its deadline '
                               '(DNS or network setup may be stuck). The import was stopped; please retry.')
        if returncode is not None:
            raise RuntimeError('The download helper stopped unexpectedly')
        now = monotonic()
        if self.disconnected is not None:
            # EOF and the OS exit status can become visible in different polls.
            # Give the specific stalled-worker status a chance to arrive without
            # blocking Blender or misreporting it as a generic disconnection.
            if now - self.disconnected > 1.:
                raise RuntimeError('The download helper disconnected')
            return None
        # Count only time spent polling responsively. Long glTF imports, modal
        # tools or debugger pauses in Blender must not expire valid replies.
        self.silent_poll_seconds += min(now - self.last_poll, .25)
        self.last_poll = now
        if self.connection is None:
            try:
                self.connection, _ = self.listener.accept()
                self.connection.setblocking(False)
            except BlockingIOError:
                return self._idle()
        deadline = monotonic() + 0.004
        # At most one response-sized frame plus a read chunk is buffered here.
        # Subsequent replies stay in the helper/socket until the caller consumes.
        while monotonic() < deadline:
            if self.authenticated and len(self.incoming) >= 17:
                kind, identity, size = struct.unpack('!cQQ', self.incoming[:17])
                if size > max(self.max_bytes, 524288):
                    raise RuntimeError('Download exceeds the response size limit')
                if len(self.incoming) >= size + 17:
                    break
            try:
                chunk = self.connection.recv(256 * 1024)
            except BlockingIOError:
                break
            except ConnectionError:
                # Windows can report a reset instead of EOF when the supervisor
                # exits the helper. Preserve its exit-code diagnosis as above.
                self.disconnected = monotonic()
                return None
            if not chunk:
                self.disconnected = monotonic()
                return None
            self.incoming.extend(chunk)
            self.silent_poll_seconds = 0.
            if not self.authenticated and len(self.incoming) >= 32:
                if bytes(self.incoming[:32]) != self.token:
                    raise RuntimeError('Invalid download helper connection')
                del self.incoming[:32]
                self.authenticated = True
                self.listener.close()
                self.listener = None
        if not self.authenticated:
            return self._idle()
        if self.outgoing:
            try:
                sent = self.connection.send(self.outgoing)
                del self.outgoing[:sent]
            except BlockingIOError:
                pass
            except ConnectionError:
                self.disconnected = monotonic()
                return None
        if len(self.incoming) < 17:
            return self._idle()
        kind, identity, size = struct.unpack('!cQQ', self.incoming[:17])
        if kind == b'H' and identity == 0 and size == 0:
            del self.incoming[:17]
            self.silent_poll_seconds = 0.
            return None
        if kind not in (b'B', b'E') or identity not in self.pending:
            raise RuntimeError('Invalid response from the download helper')
        if size > (self.max_bytes if kind == b'B' else 524288):
            raise RuntimeError('Download exceeds the response size limit')
        if len(self.incoming) < size + 17:
            return self._idle()
        content = bytes(memoryview(self.incoming)[17:17 + size])
        del self.incoming[:17 + size]
        url = self.pending.pop(identity)
        if kind == b'B':
            return identity, content, None
        detail = json.loads(content)
        if 'status' in detail:
            failure = HTTPError(redact_url(url), detail['status'], str(DownloadError(detail)),
                detail.get('headers', {}), io.BytesIO(detail.get('body', '').encode('utf-8')))
            failure.attempts = detail.get('attempts', 0)
            failure.auth_failed = detail.get('auth_failed', False)
            failure.download_code = detail.get('code', 'http_error')
        else:
            failure = DownloadError(detail)
        return identity, None, failure

    def _idle(self):
        if self.silent_poll_seconds > 15:
            raise TimeoutError('The download helper is unresponsive')
        return None

    def _check_online_access(self):
        try:
            require_online_access()
        except RuntimeError:
            self.close()
            raise

    def close(self):
        _live.discard(self)
        for connection in (self.connection, self.listener):
            if connection is not None:
                connection.close()
        self.connection = self.listener = None
        self.incoming.clear()
        self.outgoing.clear()
        self.pending.clear()
        if self.process is not None:
            if self.process.poll() is None:
                self.process.kill()
            _retired.append(self.process)
            self.process = None
            reap()
