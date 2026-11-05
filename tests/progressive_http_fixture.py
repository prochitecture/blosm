"""Isolated HTTP fixture for Blender UI tests; no fixtures are saved to disk."""
import base64
from http.server import BaseHTTPRequestHandler, HTTPServer
import json
import sys
import time

configuration = json.loads(sys.stdin.readline())
payloads = {name: base64.b64decode(value) for name, value in configuration['payloads'].items()}


class Handler(BaseHTTPRequestHandler):
    def log_message(self, *_):
        pass

    def do_GET(self):
        key = self.path.lstrip('/')
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


server = HTTPServer(('127.0.0.1', 0), Handler)
print(server.server_port, flush=True)
server.serve_forever()
