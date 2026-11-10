"""A crashed/finished Blender must not leave its HTTP or download helpers alive."""
import ctypes
import os
from pathlib import Path
import subprocess
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent))
from release_process import owned_process


@unittest.skipUnless(os.name == 'nt', 'Windows job-object contract')
class ProcessOwnershipTests(unittest.TestCase):
    def test_child_is_killed_after_parent_exits(self):
        from ctypes import wintypes
        api = ctypes.WinDLL('kernel32', use_last_error=True)
        api.OpenProcess.argtypes = [wintypes.DWORD, wintypes.BOOL, wintypes.DWORD]
        api.OpenProcess.restype = wintypes.HANDLE
        api.WaitForSingleObject.argtypes = [wintypes.HANDLE, wintypes.DWORD]
        api.WaitForSingleObject.restype = wintypes.DWORD
        api.CloseHandle.argtypes = [wintypes.HANDLE]
        script = (
            "import subprocess, sys; "
            "child = subprocess.Popen([sys.executable, '-I', '-c', 'import time; time.sleep(60)'], "
            "stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL); "
            "print(child.pid, flush=True); sys.stdin.read(1)"
        )
        handle = None
        try:
            with owned_process([sys.executable, '-I', '-c', script], stdin=subprocess.PIPE,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE) as parent:
                child_pid = int(parent.stdout.readline())
                handle = api.OpenProcess(0x100000, False, child_pid)
                self.assertTrue(handle)
                self.assertEqual(api.WaitForSingleObject(handle, 0), 258)  # WAIT_TIMEOUT: alive
                parent.stdin.write(b'\n')
                parent.stdin.flush()
                self.assertEqual(parent.wait(timeout=5), 0)
                parent.stdin.close()
                parent.stdout.close()
                parent.stderr.close()
            self.assertEqual(api.WaitForSingleObject(handle, 5000), 0, 'Orphan helper survived job close')
        finally:
            if handle:
                api.CloseHandle(handle)


if __name__ == '__main__':
    unittest.main()
