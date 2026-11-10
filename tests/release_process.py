"""Own the test subprocess tree, including helpers left behind by a Blender crash."""
from contextlib import contextmanager
import os
import signal
import subprocess


def windows_job(process):
    import ctypes as c
    from ctypes import wintypes as w

    class Basic(c.Structure):
        _fields_ = [('process_time', c.c_int64), ('job_time', c.c_int64),
                    ('flags', w.DWORD), ('minimum', c.c_size_t), ('maximum', c.c_size_t),
                    ('active', w.DWORD), ('affinity', c.c_size_t),
                    ('priority', w.DWORD), ('scheduling', w.DWORD)]

    class Io(c.Structure):
        _fields_ = [(name, c.c_uint64) for name in ('read_ops', 'write_ops', 'other_ops',
                                                  'read_bytes', 'write_bytes', 'other_bytes')]

    class Limits(c.Structure):
        _fields_ = [('basic', Basic), ('io', Io), ('process_memory', c.c_size_t),
                    ('job_memory', c.c_size_t), ('peak_process', c.c_size_t), ('peak_job', c.c_size_t)]

    api = c.WinDLL('kernel32', use_last_error=True)
    api.CreateJobObjectW.argtypes = [c.c_void_p, w.LPCWSTR]
    api.CreateJobObjectW.restype = w.HANDLE
    api.SetInformationJobObject.argtypes = [w.HANDLE, c.c_int, c.c_void_p, w.DWORD]
    api.SetInformationJobObject.restype = w.BOOL
    api.AssignProcessToJobObject.argtypes = [w.HANDLE, w.HANDLE]
    api.AssignProcessToJobObject.restype = w.BOOL
    api.CloseHandle.argtypes = [w.HANDLE]
    api.CloseHandle.restype = w.BOOL
    handle = api.CreateJobObjectW(None, None)
    if not handle:
        raise c.WinError(c.get_last_error())
    try:
        limits = Limits()
        limits.basic.flags = 0x2000  # JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE
        if not api.SetInformationJobObject(handle, 9, c.byref(limits), c.sizeof(limits)):
            raise c.WinError(c.get_last_error())
        if not api.AssignProcessToJobObject(handle, int(process._handle)):
            raise c.WinError(c.get_last_error())
    except Exception:
        api.CloseHandle(handle)
        raise
    return lambda: api.CloseHandle(handle)


@contextmanager
def owned_process(command, **kwargs):
    process = subprocess.Popen(command, start_new_session=os.name != 'nt', **kwargs)
    close = None
    try:
        if os.name == 'nt':
            close = windows_job(process)
        yield process
    finally:
        if close:
            close()
        elif os.name != 'nt':
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        elif process.poll() is None:
            process.kill()
        process.wait()
