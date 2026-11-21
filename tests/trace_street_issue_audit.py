"""Run the unchanged native audit with task tracing and source-epoch evidence.

Only the test process is instrumented. This wrapper neither provides numerical
products nor modifies source/task digests, evaluation order or acceptance.
"""
import faulthandler
import hashlib
import json
from pathlib import Path
import runpy
import sys
import time

sys.dont_write_bytecode = True
REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / 'tests'))
import check_street_progressive as support

arguments = sys.argv[sys.argv.index('--') + 1:]
OUTPUT = Path(arguments[arguments.index('--artifacts') + 1])
OUTPUT.mkdir(parents=True, exist_ok=True)
started = time.perf_counter()


def source_epoch():
    return {str(path.relative_to(REPO)).replace('\\', '/'):
            hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted(REPO.rglob('*.py'))
            if '.git' not in path.parts and '__pycache__' not in path.parts}


before = source_epoch()
original_enable = support.enable
trace_stream = (OUTPUT / 'native_task_trace.jsonl').open('w', encoding='utf-8')
stack_stream = (OUTPUT / 'native_stack_trace.log').open('w', encoding='utf-8')
faulthandler.enable(file=stack_stream)
faulthandler.dump_traceback_later(120, repeat=True, file=stack_stream)


def enable_with_trace():
    package = original_enable()
    from blosm.way.street_network.python_backend.updates import UpdateCoordinator
    original_steps = UpdateCoordinator.prepare_steps

    def traced_steps(self, network, *args, **kwargs):
        original_progress = kwargs.get('on_progress')

        def on_progress(label):
            if original_progress is not None:
                original_progress(label)
            if label.startswith('Building ') or label == 'Street geometry prepared':
                event = dict(seconds=time.perf_counter() - started,
                             network=network.uuid, label=label)
                trace_stream.write(json.dumps(event) + '\n')
                trace_stream.flush()
                print('NATIVE_TASK', json.dumps(event), flush=True)

        kwargs['on_progress'] = on_progress
        return (yield from original_steps(self, network, *args, **kwargs))

    UpdateCoordinator.prepare_steps = traced_steps
    return package


support.enable = enable_with_trace
try:
    runpy.run_path(str(REPO / 'tests/check_street_issue_audit.py'), run_name='__main__')
finally:
    faulthandler.cancel_dump_traceback_later()
    faulthandler.disable()
    trace_stream.close()
    stack_stream.close()
    after = source_epoch()
    changed = sorted(name for name in before.keys() | after.keys()
                     if before.get(name) != after.get(name))
    (OUTPUT / 'source_epoch.json').write_text(json.dumps(dict(
        changed_during_run=changed, frozen_source_epoch=not changed,
        before=before, after=after), indent=2), encoding='utf-8')
