"""Isolated test entry point for an uninterruptible DNS call and short grace."""
import os
from pathlib import Path
import socket
import sys
import threading

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from downloads import worker

mode = sys.argv.pop(1)
worker.WORKER_DEADLINE_GRACE = .35
if mode == 'dns':
    # Make the targeted getaddrinfo call independent of the test host's proxies.
    for name in list(os.environ):
        if name.lower().endswith('_proxy'):
            del os.environ[name]
    real_getaddrinfo = socket.getaddrinfo
    def getaddrinfo(host, *args, **kwargs):
        if host == 'stalled.invalid':
            threading.Event().wait(60)
        return real_getaddrinfo(host, *args, **kwargs)
    socket.getaddrinfo = getaddrinfo
elif mode != 'normal':
    raise ValueError('Unknown test mode')

worker.run()
