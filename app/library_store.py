"""Retain exact, self-contained Blender library revisions outside addon updates.

This module does not import Blender. Files are streamed in bounded chunks and
published atomically without replacing an existing revision. The caller must
certify that a library has no unpacked file dependencies before relocating it.
"""
import hashlib
import os
from pathlib import Path
import re
import tempfile


class AssetStoreError(RuntimeError):
    pass


def validate_digest(digest):
    if not isinstance(digest, str) or re.fullmatch(r'[0-9a-f]{64}', digest) is None:
        raise ValueError('A library fingerprint must be a lowercase SHA-256 digest')
    return digest


def file_digest(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def _publish(temporary, target):
    # Windows rename is atomic and refuses replacement. POSIX rename replaces
    # existing files, so use atomic no-replace hard-link creation there instead.
    if os.name == 'nt':
        os.rename(temporary, target)
    else:
        os.link(temporary, target)


class ImmutableLibraryStore:
    """Explicit durable storage, never a disposable cache or an upgrade target.

    A saved .blend may still reference any installed revision. Consequently this
    service has no eviction and no automatic in-place repair/replacement. Copy
    the store with a project when moving it; retain the digest directory layout.
    """
    def __init__(self, root):
        self.root = Path(root).resolve()

    def path(self, digest):
        return self.root / validate_digest(digest) / 'library.blend'

    def existing(self, digest):
        path = self.path(digest)
        if not path.exists():
            return None
        self._verify(path, digest)
        return path

    def install(self, source, digest):
        existing = self.existing(digest)
        if existing is not None:
            return existing
        target = self.path(digest)
        temporary = None
        try:
            target.parent.mkdir(parents=True, exist_ok=True)
            value = hashlib.sha256()
            with Path(source).open('rb') as source_stream, tempfile.NamedTemporaryFile(
                    dir=target.parent, prefix='install-', suffix='.partial', delete=False) as output:
                temporary = Path(output.name)
                for block in iter(lambda: source_stream.read(1024 * 1024), b''):
                    value.update(block)
                    output.write(block)
                output.flush()
                os.fsync(output.fileno())
            if value.hexdigest() != digest:
                raise AssetStoreError(f'Library fingerprint mismatch: {source}; expected {digest}')
            try:
                # A competing installer cannot be replaced or read half-written.
                _publish(temporary, target)
            except FileExistsError:
                self._verify(target, digest)
            return target
        except OSError as exc:
            raise AssetStoreError(f'Cannot retain library {digest} in {self.root}: {exc}') from exc
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)

    @staticmethod
    def _verify(path, digest):
        try:
            actual = file_digest(path)
        except OSError as exc:
            raise AssetStoreError(f'Cannot read retained library: {path}: {exc}') from exc
        if actual != digest:
            raise AssetStoreError(f'Retained library is damaged or was replaced: {path}; '
                                  f'expected {digest}. Restore this exact revision.')
