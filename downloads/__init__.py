"""Reusable in-memory HTTP downloads; no Blender dependencies."""
from dataclasses import dataclass, field

# A specific helper exit code preserves the reason even if IPC is backpressured
# or the caller is busy importing a previously completed result.
STALLED_WORKER_EXIT_CODE = 73


@dataclass(frozen=True)
class Request:
    url: str
    headers: dict = field(default_factory=dict)
    decode_gzip: bool = False
    verify_tls: bool = True
    timeout: float = 90.
    auth_scope: str = ''
    auth_error_reasons: tuple = ()


class DownloadError(RuntimeError):
    """A resource failure, distinct from loss of the transport/helper."""
    def __init__(self, detail):
        self.code = detail.get('code', 'download_error')
        self.attempts = detail.get('attempts', 0)
        self.auth_failed = detail.get('auth_failed', False)
        super().__init__('%s (attempts: %s)' % (detail['message'], self.attempts))
