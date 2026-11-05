"""Private transport helper. Runs in ordinary Python, never imports bpy."""
import gzip
import io
import json
import socket
import ssl
import struct
import sys
from urllib import request, error

MAX_BYTES = 256 * 1024 * 1024


def read_exact(connection, size):
    chunks = bytearray()
    while len(chunks) < size:
        chunk = connection.recv(size - len(chunks))
        if not chunk:
            raise EOFError
        chunks.extend(chunk)
    return chunks


def download(url):
    if not url.startswith(('https://', 'http://')):
        raise ValueError('Progressive downloads require an HTTP or HTTPS URL')
    req = request.Request(url, headers={'User-Agent': 'Blosm'})
    # Match the existing blocking downloader's TLS behavior.
    with request.urlopen(req, context=ssl._create_unverified_context(), timeout=30) as response:
        content = response.read(MAX_BYTES + 1)
        expected = response.headers.get('Content-Length')
        if expected is not None and len(content) < int(expected) and len(content) <= MAX_BYTES:
            raise ValueError('Incomplete HTTP response')
    if len(content) > MAX_BYTES:
        raise ValueError('Tile exceeds the 256 MiB progressive import limit')
    if content[:2] == b'\x1f\x8b':
        with gzip.GzipFile(fileobj=io.BytesIO(content)) as stream:
            content = stream.read(MAX_BYTES + 1)
        if len(content) > MAX_BYTES:
            raise ValueError('Decoded tile exceeds the 256 MiB progressive import limit')
    return content


def main():
    with socket.create_connection(('127.0.0.1', int(sys.argv[1])), timeout=15) as connection:
        connection.settimeout(None)
        connection.sendall(sys.argv[2].encode('ascii'))
        while True:
            size, = struct.unpack('!I', read_exact(connection, 4))
            if size > 65536:
                raise ValueError('Invalid helper request')
            job = json.loads(read_exact(connection, size))
            try:
                content = download(job['url'])
                kind = b'B'
            except error.HTTPError as exc:
                content = json.dumps({'status': exc.code, 'message': str(exc.reason),
                                      'body': exc.read(65536).decode('utf-8', errors='replace')}).encode()
                kind = b'E'
            except Exception as exc:
                # Avoid including the URL (which may contain API credentials).
                content = json.dumps({'message': 'Download failed: ' + type(exc).__name__}).encode()
                kind = b'E'
            connection.sendall(struct.pack('!cQ', kind, len(content)))
            connection.sendall(content)


if __name__ == '__main__':
    try:
        main()
    except (EOFError, ConnectionError):
        pass
