"""Isolated HTTP fixture for Blender UI tests; no fixtures are saved to disk."""
import base64
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json
import sys
import time
import threading

configuration = json.loads(sys.stdin.readline())
payloads = {name: base64.b64decode(value) for name, value in configuration['payloads'].items()}
guard = threading.Lock()
active = peak = 0
requests = []


class Handler(BaseHTTPRequestHandler):
    def log_message(self, *_):
        pass

    def do_GET(self):
        global active, peak
        if self.path == '/stats':
            with guard:
                content = json.dumps({'peak': peak, 'requests': requests}).encode()
            self.send_response(200)
            self.send_header('Content-Length', str(len(content)))
            self.end_headers()
            self.wfile.write(content)
            return
        with guard:
            active += 1
            peak = max(peak, active)
            requests.append(self.path)
        try:
            self.serve_payload()
        finally:
            with guard:
                active -= 1

    def serve_payload(self):
        key = self.path.split('?', 1)[0].lstrip('/')
        time.sleep(configuration.get('delays', {}).get(key, 0))
        content = payloads.get(key)
        if content is None:
            self.send_error(404)
            return
        self.send_response(200)
        self.send_header('Content-Length', str(len(content)))
        self.end_headers()
        try:
            self.wfile.write(content)
        except ConnectionError:
            pass


server = ThreadingHTTPServer(('127.0.0.1', 0), Handler)
print(server.server_port, flush=True)
server.serve_forever()
