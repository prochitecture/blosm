"""Single-request, memory-only HTTP transport in a separate Python process.

Only nonblocking sockets are used by the caller. No worker threads, dataset
files or credentials on the command line. Future transports can implement the
same submit/poll/close contract without changing Job or Blender's UI driver.
"""
import atexit
import io
import json
from pathlib import Path
import secrets
import socket
import struct
import subprocess
import sys
from time import monotonic
from urllib.error import HTTPError

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
    # Blender embeds Python; sys.executable is normally blender.exe.
    root = Path(sys.prefix)
    for candidate in (root / 'bin' / 'python.exe', root / 'bin' / 'python3',
                      root / 'bin' / ('python%d.%d' % sys.version_info[:2])):
        if candidate.is_file():
            return str(candidate)
    if Path(sys.executable).name.lower().startswith('python'):
        return sys.executable
    raise RuntimeError("Could not locate Blender's bundled Python for progressive importing")


class DownloadTransport:
    def __init__(self):
        self.listener = socket.socket()
        self.listener.bind(('127.0.0.1', 0))
        self.listener.listen(1)
        self.listener.setblocking(False)
        self.connection = None
        self.token = secrets.token_hex(16).encode('ascii')
        self.authenticated = False
        self.incoming = bytearray()
        self.outgoing = bytearray()
        self.pending = None
        self.deadline = monotonic() + 15
        self.process = None
        try:
            self.process = subprocess.Popen(
                [python_executable(), '-I', '-B', str(Path(__file__).with_name('worker.py')),
                 str(self.listener.getsockname()[1]), self.token.decode('ascii')],
                stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                creationflags=getattr(subprocess, 'CREATE_NO_WINDOW', 0))
            _live.add(self)
        except Exception:
            self.close()
            raise

    def submit(self, url):
        if self.pending is not None:
            raise RuntimeError('A download is already pending')
        data = json.dumps({'url': url}).encode('utf-8')
        if len(data) > 65536:
            raise ValueError('Download URL is too long')
        self.pending = url
        self.outgoing.extend(struct.pack('!I', len(data)) + data)
        self.deadline = monotonic() + 90

    def poll(self):
        if self.process.poll() is not None:
            raise RuntimeError('The progressive download helper stopped unexpectedly')
        if monotonic() > self.deadline:
            raise TimeoutError('The progressive download helper timed out')
        if self.connection is None:
            try:
                self.connection, _ = self.listener.accept()
                self.connection.setblocking(False)
            except BlockingIOError:
                return None
        # Bound each IPC turn as well as the adapter's work. Large replies may
        # take several turns; a partial frame is never treated as a full tile.
        deadline = monotonic() + 0.004
        while monotonic() < deadline:
            try:
                chunk = self.connection.recv(256 * 1024)
            except BlockingIOError:
                break
            if not chunk:
                raise RuntimeError('The progressive download helper disconnected')
            self.incoming.extend(chunk)
            if len(self.incoming) > MAX_BYTES + 41:
                raise RuntimeError('Download exceeds the progressive import size limit')
        if not self.authenticated:
            if len(self.incoming) < len(self.token):
                return None
            if bytes(self.incoming[:32]) != self.token:
                raise RuntimeError('Invalid progressive helper connection')
            del self.incoming[:32]
            self.authenticated = True
            self.listener.close()
            self.listener = None
        if self.outgoing:
            try:
                sent = self.connection.send(self.outgoing)
                del self.outgoing[:sent]
            except BlockingIOError:
                pass
        if len(self.incoming) < 9:
            return None
        kind, size = struct.unpack('!cQ', self.incoming[:9])
        if size > MAX_BYTES:
            raise RuntimeError('Download exceeds the progressive import size limit')
        if len(self.incoming) < size + 9:
            return None
        content = bytes(memoryview(self.incoming)[9:9 + size])
        del self.incoming[:9 + size]
        url, self.pending = self.pending, None
        if kind == b'B':
            return content, None
        if kind != b'E':
            raise RuntimeError('Invalid response from the progressive download helper')
        detail = json.loads(content)
        if 'status' in detail:
            return None, HTTPError(url, detail['status'], detail['message'], {},
                                   io.BytesIO(detail.get('body', '').encode('utf-8')))
        return None, RuntimeError(detail['message'])

    def close(self):
        _live.discard(self)
        for connection in (self.connection, self.listener):
            if connection is not None:
                connection.close()
        self.connection = self.listener = None
        self.incoming.clear()
        self.outgoing.clear()
        if self.process is not None:
            if self.process.poll() is None:
                self.process.kill()
            _retired.append(self.process)
            self.process = None
            reap()
