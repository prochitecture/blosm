"""Drive an import generator in bounded main-thread turns.

An adapter yields None to give the UI a turn, or a request understood by its
transport. Replies are sent back into the generator; failures are thrown into
it. The generator's return value is the final result. Its finally blocks must
release adapter resources when closed (cancellation or application shutdown).
"""
from time import monotonic


class Job:
    def __init__(self, steps, transport, validate=None):
        self.steps = steps
        self.transport = transport
        self.validate = validate
        self.waiting = False
        self.reply = None
        self.done = False
        self.cancelled = False
        self.result = None
        self.error = None

    def poll(self):
        """Collect I/O without advancing scene work (also useful while paused)."""
        if not self.done and self.waiting and self.reply is None:
            self.reply = self.transport.poll()

    def step(self, budget=0.008, max_steps=32):
        if self.done:
            return
        deadline = monotonic() + budget
        try:
            if self.validate is not None:
                self.validate()
            for _ in range(max_steps):
                value, error = None, None
                if self.waiting:
                    self.poll()
                    if self.reply is None:
                        return
                    value, error = self.reply
                    self.reply = None
                    self.waiting = False
                request = self.steps.throw(error) if error else self.steps.send(value)
                if request is not None:
                    self.transport.submit(request)
                    self.waiting = True
                    return
                if monotonic() >= deadline:
                    return
        except StopIteration as finished:
            self.result = finished.value
            self._finish()
        except Exception as exc:
            self.error = exc
            self._finish()

    def cancel(self):
        if not self.done:
            self.cancelled = True
            self._finish()

    def _finish(self):
        self.done = True
        self.reply = None
        try:
            self.steps.close()
        except Exception as exc:
            if self.error is None:
                self.error = exc
        finally:
            self.transport.close()
