"""Standalone local HTTP fixture for real Mapillary transport/UI checks."""
import json
from pathlib import Path
import sys
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import urlparse, parse_qs

sys.dont_write_bytecode=True
sys.path.insert(0,str(Path(__file__).resolve().parent))
from mapillary_test_support import fixture

config=json.loads(sys.stdin.readline())
tile,content,entries=fixture(config.get('count',96))
guard=threading.Lock()
active=peak=requests=0


class Handler(BaseHTTPRequestHandler):
    def log_message(self,*args): pass

    def do_GET(self):
        global active,peak,requests
        with guard:
            active+=1
            requests+=1
            peak=max(peak,active)
        try:
            query=parse_qs(urlparse(self.path).query)
            if self.path.startswith('/tiles/'):
                payload=content
                time.sleep(.08)
            elif self.path.startswith('/images'):
                ids=query['image_ids'][0].split(',')
                # Deliberately allow later batches to complete first.
                time.sleep(.12 if int(ids[0])%4 else .03)
                gate=config.get('gate')
                if gate and int(ids[0])-min(map(int,entries))>=config.get('gate_after',32):
                    deadline=time.monotonic()+60
                    while not Path(gate).exists() and time.monotonic()<deadline:
                        time.sleep(.02)
                if config.get('fail_metadata'):
                    self.send_error(404)
                    return
                payload=json.dumps({'data':[entries[i] for i in reversed(ids) if i in entries]}).encode()
            else:
                self.send_error(404)
                return
            self.send_response(200)
            self.send_header('Content-Length',str(len(payload)))
            self.end_headers()
            self.wfile.write(payload)
        except (ConnectionError,OSError):
            pass
        finally:
            with guard:
                active-=1
                if config.get('stats'):
                    path=Path(config['stats'])
                    temporary=path.with_suffix('.tmp')
                    temporary.write_text(json.dumps({'requests':requests,'peak':peak}))
                    temporary.replace(path)


with ThreadingHTTPServer(('127.0.0.1',0),Handler) as server:
    print(server.server_port,flush=True)
    server.serve_forever()
