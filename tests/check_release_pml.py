"""Run the legacy PML function tests, including both duplicate-named definitions."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import ast
import importlib
import tempfile
from types import SimpleNamespace
import unittest
from release_support import arguments, enable, REPO, PACKAGE, run

args = arguments()
enable(args)
(args.artifacts / 'temporary').mkdir(exist_ok=True)
tempfile.tempdir = str(args.artifacts / 'temporary')
PML = importlib.import_module(PACKAGE + '.pml').PML


def make_test(text, reference):
    # Test vectors stay in the checkout; only the selected runtime compiles them.
    with tempfile.TemporaryDirectory() as folder:
        path = Path(folder) / 'test.pml'
        path.write_text(text, encoding='utf-8')
        output = PML(str(path), '').getPythonCode()
    ast.parse(output)
    assert ''.join(output.split()) == ''.join(('styles = [%s]' % reference).split())


class PmlChecks(unittest.TestCase):
    def test_style_condition_uses_the_current_items_style(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'condition.pml'
            path.write_text('facade(item.front and style.count <= 1) { class: door; }', encoding='utf-8')
            output = PML(str(path), '').getPythonCode()
        namespace = {'Facade': lambda **attributes: SimpleNamespace(**attributes)}
        exec(compile(output, str(path), 'exec'), namespace)
        condition = namespace['styles'][0].condition
        self.assertTrue(condition(SimpleNamespace(front=True, styleBlock=SimpleNamespace(count=1))))
        self.assertFalse(condition(SimpleNamespace(front=True, styleBlock=SimpleNamespace(count=2))))
        self.assertFalse(condition(SimpleNamespace(front=False, styleBlock=SimpleNamespace(count=1))))


for path in sorted((REPO / 'pml/tests').glob('test_*.py')):
    tree = ast.parse(path.read_text(encoding='utf-8'), filename=str(path))
    for function in tree.body:
        if isinstance(function, ast.FunctionDef) and function.name.startswith('test_'):
            name = f'{function.name}_{path.stem}_line_{function.lineno}'
            namespace = {'makeTest': make_test}
            exec(compile(ast.Module(body=[function], type_ignores=[]), str(path), 'exec'), namespace)
            original = namespace[function.name]
            def check(self, original=original):
                original()
            setattr(PmlChecks, name, check)

run(args, PmlChecks)
