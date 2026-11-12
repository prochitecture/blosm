"""Bounded compositor work outside the interactive Blender process."""
import os
from pathlib import Path
import struct
import subprocess
import tempfile
import time
import zlib
import shutil
import bpy
from ..progressive.core import YieldToUI


def own_process_tree(process):
    """On Windows, parent death must also stop private Blender render children."""
    if os.name != 'nt':
        return None
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
    api.AssignProcessToJobObject.argtypes = [w.HANDLE, w.HANDLE]
    api.CloseHandle.argtypes = [w.HANDLE]
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


def placeholder(path, scene):
    nodes = scene.compositing_node_group.nodes if hasattr(scene, 'compositing_node_group') else scene.node_tree.nodes
    size = (1, 1)
    for name in ('empty_image', 'facade_texture', 'cladding_texture'):
        node = nodes.get(name)
        if node and getattr(node, 'image', None):
            size = tuple(node.image.size)
            break
    width, height = [max(1, int(n)) for n in size]
    def chunk(kind, data):
        return struct.pack('!I', len(data)) + kind + data + struct.pack('!I', zlib.crc32(kind+data))
    compress = zlib.compressobj()
    with open(path, 'wb') as stream:
        stream.write(b'\x89PNG\r\n\x1a\n')
        stream.write(chunk(b'IHDR', struct.pack('!2I5B', width, height, 8, 6, 0, 0, 0)))
        row = bytes(width*4+1)
        for _ in range(height):
            data = compress.compress(row)
            if data:
                stream.write(chunk(b'IDAT', data))
        stream.write(chunk(b'IDAT', compress.flush()))
        stream.write(chunk(b'IEND', b''))


class BakeQueue:
    def __init__(self, directory):
        self.directory = tempfile.TemporaryDirectory(prefix='.blosm-bakes-', dir=directory)
        self.tasks = []
        self.process = self.log = None
        self.close_tree = None
        self.pending_placeholders = []
        self.serial = 0
        profile_parent = Path(os.environ.get('BLENDER_USER_RESOURCES', self.directory.name)) / 'bake_workers'
        profile_parent.mkdir(parents=True, exist_ok=True)
        self.profiles = tempfile.TemporaryDirectory(prefix='osm-', dir=profile_parent)

    def submit(self, scene, target):
        if len(self.tasks) >= 128:
            raise RuntimeError('A building requires more than 128 uncached textures')
        folder = Path(self.directory.name) / str(self.serial)
        self.serial += 1
        folder.mkdir()
        snapshot = folder / 'input.blend'
        bpy.data.libraries.write(str(snapshot), {scene}, path_remap='ABSOLUTE')
        placeholder(target, scene)
        self.pending_placeholders.append(target)
        self.tasks.append((folder, target))

    def remove_placeholders(self):
        for path in self.pending_placeholders:
            if os.path.isfile(path):
                os.unlink(path)
        self.pending_placeholders.clear()

    def steps(self):
        self.remove_placeholders()
        try:
            for folder, target in self.tasks:
                output = folder / 'output'
                output.mkdir()
                profile = Path(self.profiles.name) / folder.name
                env = os.environ.copy()
                for kind in ('RESOURCES', 'CONFIG', 'SCRIPTS', 'EXTENSIONS', 'DATAFILES'):
                    location = profile / kind.lower()
                    location.mkdir(parents=True)
                    env['BLENDER_USER_'+kind] = str(location)
                self.log = (folder / 'render.log').open('wb')
                self.process = subprocess.Popen([bpy.app.binary_path, '--background',
                    '--factory-startup', '--disable-autoexec', '--python-exit-code', '1',
                    '--python', str(Path(__file__).with_name('bake_worker.py')), '--',
                    str(folder / 'input.blend'), str(output)], env=env,
                    stdin=subprocess.DEVNULL, stdout=self.log, stderr=subprocess.STDOUT,
                    creationflags=getattr(subprocess, 'CREATE_NO_WINDOW', 0))
                self.close_tree = own_process_tree(self.process)
                deadline = time.monotonic()+300
                while self.process.poll() is None:
                    if time.monotonic() >= deadline:
                        raise TimeoutError('Texture render exceeded five minutes')
                    yield YieldToUI()
                self.log.close()
                self.log = None
                if self.process.returncode:
                    detail = (folder/'render.log').read_text(errors='replace')[-1500:]
                    raise RuntimeError('Texture rendering failed: '+detail)
                files = list(output.iterdir())
                if len(files) != 1 or not files[0].is_file():
                    raise RuntimeError('Texture rendering did not produce exactly one file')
                os.replace(files[0], target)
                for image in bpy.data.images:
                    if image.source == 'FILE' and os.path.realpath(bpy.path.abspath(image.filepath)) == os.path.realpath(target):
                        image.reload()
                self.process = None
                if self.close_tree:
                    self.close_tree()
                    self.close_tree = None
                if folder.resolve().parent != Path(self.directory.name).resolve():
                    raise RuntimeError('Invalid private bake directory')
                shutil.rmtree(folder)
            self.tasks.clear()
        finally:
            if self.close_tree:
                self.close_tree()
                self.close_tree = None
            if self.process and self.process.poll() is None:
                self.process.kill()
                self.process.wait(timeout=2)
            self.process = None
            if self.log:
                self.log.close()
                self.log = None

    def close(self):
        self.remove_placeholders()
        if self.close_tree:
            self.close_tree()
            self.close_tree = None
        if self.process and self.process.poll() is None:
            self.process.kill()
            self.process.wait(timeout=2)
        if self.log:
            self.log.close()
        self.directory.cleanup()
        self.profiles.cleanup()
