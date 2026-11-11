"""Optional, process-local execution diagnostics for street benchmarks.

Read Windows processor topology, current main-thread CPU and power counters.
ProcessorState only reads information. ProcessQoS is an explicit test-only,
current-process override which restores the original policy on close. No
priority, affinity, power-plan or persistent environment changes are made here.
Coarse numerical wrappers count calls and time their unchanged implementation;
short calibration kernels distinguish execution rate from construction choices.
"""
from collections import defaultdict
from contextlib import ExitStack
import ctypes
from ctypes import wintypes
import json
import os
import struct
import sys
import time
from unittest.mock import patch


class ProcessorState:
    def __init__(self):
        self.topology = []
        self.error = None
        if os.name != 'nt':
            self.error = 'Windows processor diagnostics unavailable'
            return
        self.kernel = ctypes.WinDLL('kernel32', use_last_error=True)
        self.kernel.GetCurrentProcess.restype = wintypes.HANDLE
        self.kernel.GetCurrentThread.restype = wintypes.HANDLE
        self.kernel.GetPriorityClass.argtypes = [wintypes.HANDLE]
        self.kernel.GetPriorityClass.restype = wintypes.DWORD
        self.kernel.GetCurrentProcessorNumber.restype = wintypes.DWORD
        self.kernel.GetThreadPriority.argtypes = [wintypes.HANDLE]
        self.kernel.GetThreadPriority.restype = ctypes.c_int
        self.kernel.GetSystemCpuSetInformation.argtypes = [ctypes.c_void_p, wintypes.ULONG,
            ctypes.POINTER(wintypes.ULONG), wintypes.HANDLE, wintypes.ULONG]
        size = wintypes.ULONG()
        self.kernel.GetSystemCpuSetInformation(None, 0, ctypes.byref(size), None, 0)
        buffer = ctypes.create_string_buffer(size.value)
        if self.kernel.GetSystemCpuSetInformation(buffer, size, ctypes.byref(size), None, 0):
            offset = 0
            while offset < size.value:
                length, kind = struct.unpack_from('<II', buffer, offset)
                if length < 8 or offset + length > size.value:
                    raise ValueError('Invalid CPU-set information size')
                if kind == 0 and length >= 24:
                    identity, group, logical, core, cache, numa, efficiency, flags = struct.unpack_from(
                        '<IH6B', buffer, offset + 8)
                    self.topology.append(dict(id=identity, group=group, cpu=logical, core=core,
                        cache=cache, numa=numa, efficiency_class=efficiency, flags=flags))
                offset += length
        else:
            self.error = 'GetSystemCpuSetInformation: ' + str(ctypes.get_last_error())

    def current(self):
        if os.name != 'nt':
            return {}
        cpu = self.kernel.GetCurrentProcessorNumber()
        info = next((t for t in self.topology if t['group'] == 0 and t['cpu'] == cpu), {})
        return dict(cpu=cpu, efficiency_class=info.get('efficiency_class'),
            process_priority_class=self.kernel.GetPriorityClass(self.kernel.GetCurrentProcess()),
            thread_priority=self.kernel.GetThreadPriority(self.kernel.GetCurrentThread()))

    def power(self):
        if os.name != 'nt':
            return {}
        class Power(ctypes.Structure):
            _fields_ = [(n, wintypes.ULONG) for n in
                ('number', 'maximum_mhz', 'current_mhz', 'limit_mhz', 'max_idle', 'current_idle')]
        values = (Power * (os.cpu_count() or 1))()
        api = ctypes.WinDLL('powrprof').CallNtPowerInformation
        api.argtypes = [ctypes.c_int, ctypes.c_void_p, wintypes.ULONG, ctypes.c_void_p, wintypes.ULONG]
        api.restype = wintypes.ULONG
        status = api(11, None, 0, values, ctypes.sizeof(values))
        result = dict(status=status, processors=[{n: getattr(v, n) for n, _ in Power._fields_} for v in values]
            if not status else [])
        class Throttling(ctypes.Structure):
            _fields_ = [(n, wintypes.ULONG) for n in ('version', 'control_mask', 'state_mask')]
        for function, handle, kind, label in (
                ('GetProcessInformation', self.kernel.GetCurrentProcess(), 4, 'process_power_policy'),
                ('GetThreadInformation', self.kernel.GetCurrentThread(), 3, 'thread_power_policy')):
            value = Throttling(1, 0, 0)
            query = getattr(self.kernel, function)
            query.argtypes = [wintypes.HANDLE, ctypes.c_int, ctypes.c_void_p, wintypes.DWORD]
            success = query(handle, kind, ctypes.byref(value), ctypes.sizeof(value))
            result[label] = dict(success=bool(success), error=0 if success else ctypes.get_last_error(),
                                 **{n: getattr(value, n) for n, _ in Throttling._fields_})
        return result


class ProcessQoS:
    """Explicit scheduling experiment, restricted to this test process."""
    def __init__(self, mode):
        if mode not in ('high', 'eco') or os.name != 'nt':
            raise ValueError('A Windows test process is required for high/eco QoS')
        self.state = ProcessorState()
        self.original = self.state.power()['process_power_policy']
        if not self.original['success']:
            raise OSError(self.original['error'], 'Cannot read original test-process QoS')
        self.report = dict(requested=mode, before=self.original)
        self.closed = False
        self._set(self.original['control_mask'] | 1,
                  (self.original['state_mask'] & ~1) | (mode == 'eco'))
        try:
            self.report['during'] = self.state.power()['process_power_policy']
            if not self.report['during']['success']:
                raise OSError(self.report['during']['error'], 'Cannot verify test-process QoS')
        except BaseException:
            self.close()
            raise

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
        return False

    def _set(self, control, state):
        class Policy(ctypes.Structure):
            _fields_ = [(n, wintypes.ULONG) for n in ('version', 'control', 'state')]
        value = Policy(1, control, state)
        api = self.state.kernel.SetProcessInformation
        api.argtypes = [wintypes.HANDLE, ctypes.c_int, ctypes.c_void_p, wintypes.DWORD]
        if not api(self.state.kernel.GetCurrentProcess(), 4, ctypes.byref(value), ctypes.sizeof(value)):
            raise ctypes.WinError(ctypes.get_last_error())

    def close(self):
        if self.closed:
            return
        self._set(self.original['control_mask'], self.original['state_mask'])
        self.report['after_restore'] = self.state.power()['process_power_policy']
        assert self.report['after_restore'] == self.original, self.report
        self.closed = True


class ExecutionDiagnostics:
    def __init__(self, output):
        self.output = output
        self.processor = ProcessorState()
        self.calls = defaultdict(lambda: dict(calls=0, wall_seconds=0., thread_cpu_seconds=0.,
                                              shapes=defaultdict(int), cpu_pairs=defaultdict(int)))
        self.calibrations = []
        self.stack = ExitStack()

    def calibrate(self, label):
        import numpy as np
        before = self.processor.current()
        start = time.perf_counter()
        value = 0
        for i in range(200000):
            value += (i * 17) % 1009
        python_seconds = time.perf_counter() - start
        vector = np.linspace(.1, 5., 16384).reshape(-1, 2)
        start = time.perf_counter()
        for _ in range(100):
            value = float(np.linalg.norm(vector, axis=1).sum())
        numpy_seconds = time.perf_counter() - start
        self.calibrations.append(dict(label=label, python_seconds=python_seconds,
            numpy_seconds=numpy_seconds, check=value, before=before, after=self.processor.current(),
            power=self.processor.power(), trace=repr(sys.gettrace()), profile=repr(sys.getprofile())))
        self.save()

    def wrapper(self, original, name):
        def call(*args, **kwargs):
            before = self.processor.current()
            shape = next((str(a.shape) for a in args if hasattr(a, 'shape')), '')
            start, cpu = time.perf_counter(), time.thread_time()
            try:
                return original(*args, **kwargs)
            finally:
                seconds, cpu_seconds = time.perf_counter() - start, time.thread_time() - cpu
                after = self.processor.current()
                entry = self.calls[name]
                entry['calls'] += 1
                entry['wall_seconds'] += seconds
                entry['thread_cpu_seconds'] += cpu_seconds
                entry['shapes'][shape] += 1
                pair = (before.get('cpu'), before.get('efficiency_class'), after.get('cpu'), after.get('efficiency_class'))
                entry['cpu_pairs'][str(pair)] += 1
        return call

    def install(self, package):
        import importlib
        prefix = package + '.way.street_network.python_backend.'
        targets = [('quad_charts', 'FixedRimQuality', 'improve'),
            ('quad_workspace', 'QuadEnergy', 'evaluate'), ('quad_workspace', 'QuadGradient', 'evaluate'),
            ('surface_fields', 'BoundaryJetField', '__call__'),
            ('decomposition', 'VisibilitySplitter', 'split'),
            ('partition', 'CompoundQuadChart', 'build')]
        for module, cls, method in targets:
            owner = getattr(importlib.import_module(prefix + module), cls)
            name = cls + '.' + method
            self.stack.enter_context(patch.object(owner, method, self.wrapper(getattr(owner, method), name)))
        charts = importlib.import_module(prefix + 'quad_charts')
        original = charts.visibility_center
        wrapped = self.wrapper(original, 'visibility_center')
        # Direct imports already bound aliases; patch those exact aliases too.
        for module in list(sys.modules.values()):
            if (getattr(module, '__name__', '').startswith(prefix)
                    and getattr(module, 'visibility_center', None) is original):
                self.stack.enter_context(patch.object(module, 'visibility_center', wrapped))

    def save(self):
        self.output.mkdir(parents=True, exist_ok=True)
        report = dict(topology=self.processor.topology, topology_error=self.processor.error,
            calibrations=self.calibrations, calls=dict(self.calls),
            qualification='Coarse wrappers and untimed calibration kernels; not an uninstrumented baseline')
        (self.output / 'execution_diagnostics.json').write_text(json.dumps(report, indent=2), encoding='utf8')

    def close(self):
        self.stack.close()
        self.save()


if __name__ == '__main__':
    state = ProcessorState()
    print(json.dumps(dict(topology=state.topology, error=state.error, current=state.current(), power=state.power()), indent=2))
