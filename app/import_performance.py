"""Best-effort, process-local performance preference for Import work.

Only the execution-speed QoS bit is requested. Explicit EcoQoS is respected,
and the complete previous policy is restored on leaving each work scope.
The module has no Blender dependency and loads the Windows API lazily.
"""
import ctypes
from ctypes import wintypes
import logging
import os


_EXECUTION_SPEED = 1
_logger = logging.getLogger(__name__)


class _PowerPolicy(ctypes.Structure):
    _fields_ = [('version', wintypes.ULONG), ('control', wintypes.ULONG),
                ('state', wintypes.ULONG)]


class _WindowsProcessPolicy:
    def __init__(self):
        api = ctypes.WinDLL('kernel32', use_last_error=True)
        api.GetCurrentProcess.restype = wintypes.HANDLE
        api.GetProcessInformation.argtypes = [wintypes.HANDLE, ctypes.c_int,
                                              ctypes.c_void_p, wintypes.DWORD]
        api.GetProcessInformation.restype = wintypes.BOOL
        api.SetProcessInformation.argtypes = [wintypes.HANDLE, ctypes.c_int,
                                              ctypes.c_void_p, wintypes.DWORD]
        api.SetProcessInformation.restype = wintypes.BOOL
        self.api = api
        self.process = api.GetCurrentProcess()

    def read(self):
        value = _PowerPolicy(1, 0, 0)
        if not self.api.GetProcessInformation(self.process, 4, ctypes.byref(value), ctypes.sizeof(value)):
            raise ctypes.WinError(ctypes.get_last_error())
        return value.control, value.state

    def write(self, policy):
        value = _PowerPolicy(1, *policy)
        if not self.api.SetProcessInformation(self.process, 4, ctypes.byref(value), ctypes.sizeof(value)):
            raise ctypes.WinError(ctypes.get_last_error())


class ImportPerformance:
    """A reversible HighQoS request; failures never replace an Import result."""
    def __init__(self, enabled=True):
        self.enabled = enabled
        self.applied = False
        self.status = 'disabled'
        self._policy = None
        self._original = None
        self._restore_needed = False

    def __enter__(self):
        if not self.enabled:
            return self
        if os.name != 'nt':
            self.status = 'unsupported'
            return self
        try:
            self._policy = _WindowsProcessPolicy()
            self._original = self._policy.read()
            control, state = self._original
            if control & state & _EXECUTION_SPEED:
                self.status = 'power_saving'
                return self
            requested = control | _EXECUTION_SPEED, state & ~_EXECUTION_SPEED
            if requested == self._original:
                self.status = 'already_high'
                return self
            # A failing set/read-back may still have changed the policy.
            self._restore_needed = True
            self._policy.write(requested)
            if self._policy.read() != requested:
                raise OSError('Windows did not confirm the Import performance request')
            self.applied = True
            self.status = 'active'
        except BaseException as exc:
            self.status = 'unavailable'
            self.close()
            if not isinstance(exc, (OSError, AttributeError)):
                raise
            _logger.debug('Import performance request unavailable: %s', exc)
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
        return False

    def close(self):
        if not self._restore_needed:
            return
        try:
            self._policy.write(self._original)
            if self._policy.read() != self._original:
                raise OSError('Windows did not confirm restoration of the Import power policy')
        except (OSError, AttributeError) as exc:
            self.status = 'restore_failed'
            _logger.warning('Could not restore the previous Import power policy: %s', exc)
        else:
            self._restore_needed = False
            if self.applied:
                self.status = 'restored'


def import_performance(context, addon_name):
    """Resolve the shared addon preference for synchronous or timer work."""
    entry = context.preferences.addons.get(addon_name)
    enabled = getattr(entry.preferences, 'preferPerformanceDuringImport', True) if entry else True
    return ImportPerformance(enabled)
