"""Bounded main-thread import scheduling over a reusable download transport.

Yield a URL/Request to await bytes, Download to prefetch with a main-thread
consumer, None to yield UI time, or Drain before finalizing deferred work.
"""
from collections import deque
from collections.abc import Iterator
from dataclasses import dataclass
from time import monotonic
try:
    from ..downloads import Request
except ImportError:  # Standalone transport tests/tools.
    from downloads import Request


@dataclass(frozen=True)
class Download:
    request: object
    consume: object  # Callable(data, error), optionally returning a checkpoint iterator.


class Drain:
    """Wait for all deferred consumers before resuming the generator."""


class YieldToUI:
    """End this turn immediately, including after a local multi-object commit."""


class LazyTransport:
    """Local work needs no process; allocate a transport on the first request."""
    def __init__(self, factory):
        self.factory = factory
        self.instance = None

    @property
    def max_requests(self):
        return self.instance.max_requests if self.instance else 1

    def submit(self, request):
        if self.instance is None:
            self.instance = self.factory()
        return self.instance.submit(request)

    def poll(self):
        return self.instance.poll() if self.instance else None

    def close(self):
        if self.instance:
            self.instance.close()


class Job:
    def __init__(self, steps, transport, validate=None):
        self.steps = steps
        self.transport = transport
        self.validate = validate
        self.pending = {}
        self.replies = deque()
        self.request = None
        self.awaited = None
        self.draining = False
        self.exhausted = False
        self.done = False
        self.cancelled = False
        self.result = None
        self.error = None
        self.consumer = None
        self.value = self.failure = None

    @property
    def waiting(self):
        return bool(self.pending or self.request is not None)

    def poll(self):
        """Receive while paused, without consuming results or discovering work."""
        if not self.done and len(self.replies) < len(self.pending):
            reply = self.transport.poll()
            if reply is not None:
                self.replies.append(reply)

    def step(self, budget=0.008, max_steps=32):
        if self.done:
            return
        deadline = monotonic() + budget
        try:
            if self.validate is not None:
                self.validate()
            for _ in range(max_steps):
                self.poll()
                # Parsing/committing a response may itself be cooperative. Keep
                # one active consumer, bounded by the same turn budget, before
                # admitting more decoded work. poll() never advances this iterator.
                if self.consumer is not None:
                    try:
                        next(self.consumer)
                    except StopIteration:
                        self.consumer = None
                    if monotonic() >= deadline:
                        return
                    continue
                # Refill a slot freed by the previous consumer before importing
                # the next ready tile. Otherwise fast downloads and slow imports
                # would accidentally form batches with idle workers between them.
                can_fill = len(self.pending) < self.transport.max_requests and (
                    self.request is not None or
                    (self.awaited is None and not self.draining and not self.exhausted))
                if self.replies and not can_fill:
                    identity, data, error = self.replies.popleft()
                    command = self.pending.pop(identity)
                    if isinstance(command, Download):
                        result = command.consume(data, error)
                        if isinstance(result, Iterator):
                            self.consumer = result
                        # Let the viewport display each committed result before
                        # another import or finalization, even for tiny tiles.
                        return
                    else:
                        self.awaited = None
                        self.value, self.failure = data, error
                    # Avoid retaining the last consumed tile's bytes while waiting.
                    del data, error, command
                    if monotonic() >= deadline:
                        return
                if self.request is not None:
                    if len(self.pending) >= self.transport.max_requests:
                        return
                    command = self.request
                    identity = self.transport.submit(command.request if isinstance(command, Download) else command)
                    self.pending[identity] = command
                    self.request = None
                    if not isinstance(command, Download):
                        self.awaited = identity
                    # Flush submissions now so downloads overlap the next import.
                    self.poll()
                if self.awaited is not None:
                    return
                if self.draining or self.exhausted:
                    if self.pending:
                        return
                    self.draining = False
                    if self.exhausted:
                        self._finish()
                        return
                value, failure = self.value, self.failure
                self.value = self.failure = None
                try:
                    request = self.steps.throw(failure) if failure else self.steps.send(value)
                except StopIteration as finished:
                    self.result = finished.value
                    self.exhausted = True
                    request = None
                del value, failure
                if isinstance(request, YieldToUI):
                    return
                if isinstance(request, Drain):
                    self.draining = True
                elif request is not None:
                    if not isinstance(request, (str, Request, Download)):
                        raise TypeError('Unsupported progressive command: %s' % type(request).__name__)
                    self.request = request
                if monotonic() >= deadline:
                    return
        except Exception as exc:
            self.error = exc
            self._finish()

    def cancel(self):
        if not self.done:
            self.cancelled = True
            self._finish()

    def _finish(self):
        self.done = True
        self.pending.clear()
        self.replies.clear()
        self.request = self.value = self.failure = None
        try:
            try:
                if self.consumer is not None:
                    consumer, self.consumer = self.consumer, None
                    close = getattr(consumer, 'close', None)
                    if close is not None:
                        close()
            finally:
                self.steps.close()
        except Exception as exc:
            if self.error is None:
                self.error = exc
        finally:
            self.transport.close()
