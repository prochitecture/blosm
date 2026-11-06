"""HTTP worker pool in an ordinary Python process; never imports bpy."""
from concurrent.futures import ThreadPoolExecutor
from email.utils import parsedate_to_datetime
import gzip
import io
import json
import random
import select
import socket
import ssl
import struct
import sys
import time
from urllib import request, error


def read_exact(connection, size):
    chunks = bytearray()
    while len(chunks) < size:
        chunk = connection.recv(size - len(chunks))
        if not chunk:
            raise EOFError
        chunks.extend(chunk)
    return chunks


def retry_delay(exc, attempt):
    value = exc.headers.get('Retry-After') if isinstance(exc, error.HTTPError) else None
    if value:
        try:
            return max(0., float(value))
        except ValueError:
            try:
                return max(0., parsedate_to_datetime(value).timestamp() - time.time())
            except (ValueError, TypeError, OverflowError):
                pass
    return .25 * 2 ** attempt + random.uniform(0, .1)


def download(job, max_bytes):
    url = job['url']
    if not url.startswith(('https://', 'http://')):
        raise ValueError('Downloads require HTTP or HTTPS')
    headers = {'User-Agent': 'Blosm', **job['headers']}
    ctx = ssl.create_default_context() if job['verify_tls'] else ssl._create_unverified_context()
    deadline = time.monotonic() + 80
    for attempt in range(3):
        try:
            with request.urlopen(request.Request(url, headers=headers), context=ctx, timeout=30) as response:
                content = response.read(max_bytes + 1)
                expected = response.headers.get('Content-Length')
                if expected is not None and len(content) < int(expected) and len(content) <= max_bytes:
                    raise ValueError('Incomplete HTTP response')
            break
        except (error.URLError, OSError) as exc:
            retryable = not isinstance(exc, error.HTTPError) or exc.code in (429, 500, 502, 503, 504)
            delay = retry_delay(exc, attempt)
            if not retryable or attempt == 2 or delay > 30 or time.monotonic() + delay >= deadline:
                raise
            if isinstance(exc, error.HTTPError):
                exc.close()
            time.sleep(delay)
    if len(content) > max_bytes:
        raise ValueError('Response exceeds the download size limit')
    if job['decode_gzip'] and content[:2] == b'\x1f\x8b':
        with gzip.GzipFile(fileobj=io.BytesIO(content)) as stream:
            content = stream.read(max_bytes + 1)
        if len(content) > max_bytes:
            raise ValueError('Decoded response exceeds the download size limit')
    return content


def result(job, max_bytes):
    try:
        return b'B', download(job, max_bytes)
    except error.HTTPError as exc:
        try:
            detail = {'status': exc.code, 'message': 'HTTP %s' % exc.code,
                      'headers': {'Retry-After': exc.headers.get('Retry-After', '')},
                      'body': exc.read(65536).decode('utf-8', errors='replace')}
        finally:
            exc.close()
    except Exception as exc:
        # Do not include request URLs/headers or exception text with credentials.
        detail = {'message': 'Download failed: ' + type(exc).__name__}
    return b'E', json.dumps(detail).encode()


def main():
    workers, max_bytes = int(sys.argv[3]), int(sys.argv[4])
    with socket.create_connection(('127.0.0.1', int(sys.argv[1])), timeout=15) as connection:
        connection.settimeout(None)
        connection.sendall(sys.argv[2].encode('ascii'))
        with ThreadPoolExecutor(max_workers=workers) as pool:
            pending = {}
            while True:
                if select.select([connection], [], [], .005)[0]:
                    size, = struct.unpack('!I', read_exact(connection, 4))
                    if size > 65536 or len(pending) >= workers:
                        raise ValueError('Invalid helper request')
                    job = json.loads(read_exact(connection, size))
                    pending[pool.submit(result, job, max_bytes)] = job['id']
                for future in list(pending):
                    if future.done():
                        identity = pending.pop(future)
                        kind, content = future.result()
                        connection.sendall(struct.pack('!cQQ', kind, identity, len(content)))
                        connection.sendall(content)
                        del content, future


if __name__ == '__main__':
    try:
        main()
    except (EOFError, ConnectionError):
        pass
