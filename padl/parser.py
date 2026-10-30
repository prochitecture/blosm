"""PADL parser and structural validation for the initial building compiler profile."""
import hashlib
import json
from pathlib import Path
from types import MappingProxyType

from .jsonc import JsonReader
from .model import Choice, Condition, LevelSpec, Operation, Program, Value, VolumeSpec


SCOPES = ('global', 'footprint', 'footprint_edge')
ARITIES = {'add': (2, 2), 'subtract': (2, 2), 'multiply': (2, 2), 'divide': (2, 2),
           'min': (2, None), 'max': (2, None), 'negate': (1, 1), 'abs': (1, 1),
           'floor': (1, 1), 'ceil': (1, 1)}


class PadlParser:
    """Register additional operation parsers without changing the document reader.

    A handler receives ``(parser, params_dict, params_path)`` and returns its
    backend-independent operation payload. Backend support is registered separately.
    """
    def __init__(self, operation_parsers=None):
        self.operation_parsers = {'choose_one': PadlParser.choices,
                                  'generate_volume': PadlParser.volume}
        self.operation_parsers.update(operation_parsers or {})

    def load(self, path):
        path = Path(path)
        return self.loads(path.read_text(encoding='utf-8-sig'), filename=str(path))

    def parse(self, document, filename='<mapping>'):
        return self.loads(json.dumps(document, allow_nan=False), filename=filename)

    def loads(self, text, filename='<string>'):
        document, self.source = JsonReader(text, filename).read()
        self.object(document, (), {'schema_version', 'metadata', 'context', 'definitions',
                                   'operations', 'libraries', 'footprints'})
        if document.get('schema_version') != '0.1':
            raise self.source.error(('schema_version',), 'Supported PADL schema_version is "0.1"', 'unsupported_schema')
        metadata = document.get('metadata', {})
        self.object(metadata, ('metadata',))
        if 'name' in metadata:
            self.name(metadata['name'], ('metadata', 'name'))
        for key in ('libraries', 'footprints'):
            if document.get(key):
                raise self.source.error((key,), f'Authored {key} are not supported by this external-footprint profile', 'unsupported_feature')
        context = document.get('context', {})
        self.object(context, ('context',), {'seed', 'units'})
        if context.get('units', 'm') not in ('m', 'meters', 'metres'):
            raise self.source.error(('context', 'units'), 'This backend profile requires metre units', 'unsupported_feature')
        if 'seed' in context and (type(context['seed']) is not int or not -(2**31) <= context['seed'] < 2**31):
            raise self.source.error(('context', 'seed'), 'Seed must be a signed 32-bit integer')
        definitions = document.get('definitions', {})
        self.object(definitions, ('definitions',), {'parameters', 'conditions'})
        params = definitions.get('parameters', {})
        conditions = definitions.get('conditions', {})
        self.object(params, ('definitions', 'parameters'))
        self.object(conditions, ('definitions', 'conditions'))
        parsed_params = {k: self.value(v, ('definitions', 'parameters', k)) for k, v in params.items()}
        parsed_conditions = {k: self.condition(v, ('definitions', 'conditions', k)) for k, v in conditions.items()}
        operations = self.operations(document.get('operations'), ('operations',))
        fingerprint = hashlib.sha256(json.dumps(document, sort_keys=True, separators=(',', ':'),
                                                ensure_ascii=False).encode('utf-8')).hexdigest()
        program = Program(MappingProxyType(parsed_params), MappingProxyType(parsed_conditions), operations,
                          self.source, fingerprint, document['schema_version'], MappingProxyType(context), document)
        self.check_conditions(program)
        return program

    def object(self, obj, path, allowed=None, required=()):
        if not isinstance(obj, dict):
            raise self.source.error(path, 'Expected an object')
        if allowed is not None:
            unknown = obj.keys() - allowed
            if unknown:
                key = sorted(unknown)[0]
                raise self.source.error(path + (key,), f'Unsupported field {key!r}', 'unsupported_feature')
        missing = set(required) - obj.keys()
        if missing:
            raise self.source.error(path, f'Missing required field {sorted(missing)[0]!r}')

    def array(self, obj, path, minimum=0):
        if not isinstance(obj, list) or len(obj) < minimum:
            raise self.source.error(path, f'Expected an array containing at least {minimum} item(s)')

    def name(self, value, path):
        if not isinstance(value, str) or not value:
            raise self.source.error(path, 'Expected a nonempty string')
        return value

    def value(self, obj, path):
        if isinstance(obj, (str, bool, int, float)):
            return Value('literal', obj, path)
        self.object(obj, path, required=('type',))
        kind = obj['type']
        if kind == 'constant':
            self.object(obj, path, {'type', 'value'}, ('value',))
            if not isinstance(obj['value'], (str, bool, int, float)):
                raise self.source.error(path + ('value',), 'Only scalar constants are supported', 'unsupported_feature')
            return Value('literal', obj['value'], path)
        if kind == 'parameter':
            self.object(obj, path, {'type', 'name', 'scope'}, ('name',))
            name = self.name(obj['name'], path + ('name',))
            scope = obj.get('scope', 'global')
            if scope not in SCOPES:
                raise self.source.error(path + ('scope',), f'Unknown scope {scope!r}')
            return Value(kind, (scope, name), path)
        if kind == 'expression':
            self.object(obj, path, {'type', 'op', 'args'}, ('op', 'args'))
            op = obj['op']
            if not isinstance(op, str) or op not in ARITIES:
                raise self.source.error(path + ('op',), f'Unknown arithmetic operator {op!r}')
            low, high = ARITIES[op]
            self.array(obj['args'], path + ('args',), low)
            if high is not None and len(obj['args']) > high:
                raise self.source.error(path + ('args',), f'{op} expects exactly {high} argument(s)')
            args = tuple(self.value(v, path + ('args', i)) for i, v in enumerate(obj['args']))
            return Value(kind, (op, args), path)
        if kind == 'uniform':
            self.object(obj, path, {'type', 'min', 'max'}, ('min', 'max'))
            return Value(kind, tuple(self.value(obj[k], path + (k,)) for k in ('min', 'max')), path)
        if kind in ('choice', 'weighted_choice'):
            self.object(obj, path, {'type', 'values'}, ('values',))
            self.array(obj['values'], path + ('values',), 1)
            items = []
            for i, item in enumerate(obj['values']):
                p = path + ('values', i)
                if kind == 'weighted_choice':
                    self.object(item, p, {'value', 'weight'}, ('value', 'weight'))
                    items.append((self.value(item['value'], p + ('value',)), self.value(item['weight'], p + ('weight',))))
                else:
                    items.append((self.value(item, p), Value('literal', 1, p)))
            return Value('weighted_choice', tuple(items), path)
        raise self.source.error(path + ('type',), f'Unsupported ValueSpec {kind!r}', 'unsupported_feature')

    def condition(self, obj, path):
        self.object(obj, path, required=('type',))
        kind = obj['type']
        if kind == 'ref':
            self.object(obj, path, {'type', 'name'}, ('name',))
            return Condition(kind, self.name(obj['name'], path + ('name',)), path)
        if kind in ('and', 'or'):
            self.object(obj, path, {'type', 'conditions'}, ('conditions',))
            self.array(obj['conditions'], path + ('conditions',), 1)
            return Condition(kind, tuple(self.condition(c, path + ('conditions', i)) for i, c in enumerate(obj['conditions'])), path)
        if kind == 'not':
            self.object(obj, path, {'type', 'condition'}, ('condition',))
            return Condition(kind, self.condition(obj['condition'], path + ('condition',)), path)
        if kind == 'compare':
            op = obj.get('op')
            if op not in ('eq', 'ne', 'lt', 'lte', 'gt', 'gte', 'between'):
                raise self.source.error(path + ('op',), f'Unknown comparison {op!r}')
            fields = ('left', 'min', 'max') if op == 'between' else ('left', 'right')
            self.object(obj, path, {'type', 'op', *fields}, ('op', *fields))
            return Condition(kind, (op, tuple(self.value(obj[k], path + (k,)) for k in fields)), path)
        raise self.source.error(path, f'Unsupported condition {kind!r}', 'unsupported_feature')

    def operations(self, items, path):
        self.array(items, path)
        result = []
        for i, obj in enumerate(items):
            p = path + (i,)
            self.object(obj, p, {'type', 'name', 'condition', 'params'}, ('type', 'params'))
            kind = self.name(obj['type'], p + ('type',))
            if kind not in self.operation_parsers:
                raise self.source.error(p + ('type',), f'Unsupported operation {kind!r}', 'unsupported_operation')
            cond = self.condition(obj['condition'], p + ('condition',)) if 'condition' in obj else None
            payload = self.operation_parsers[kind](self, obj['params'], p + ('params',))
            name = self.name(obj.get('name', kind), p + ('name',))
            result.append(Operation(kind, name, payload, cond, p))
        return tuple(result)

    def choices(self, obj, path):
        self.object(obj, path, {'mode', 'choices'}, ('choices',))
        if obj.get('mode', 'first_match') != 'first_match':
            raise self.source.error(path + ('mode',), 'Only first_match choices are supported', 'unsupported_feature')
        self.array(obj['choices'], path + ('choices',), 1)
        choices = []
        for i, choice in enumerate(obj['choices']):
            p = path + ('choices', i)
            self.object(choice, p, {'name', 'condition', 'operations'}, ('operations',))
            condition = self.condition(choice['condition'], p + ('condition',)) if 'condition' in choice else None
            if condition is None and i != len(obj['choices']) - 1:
                raise self.source.error(p, 'An unconditional first_match choice must be last')
            choices.append(Choice(self.name(choice.get('name', str(i)), p + ('name',)), condition,
                                  self.operations(choice['operations'], p + ('operations',)), p))
        return tuple(choices)

    def policy(self, obj, path):
        policy = obj.get('height_type', 'fixed')
        if policy not in ('fixed', 'proportional'):
            raise self.source.error(path + ('height_type',), 'Expected fixed or proportional')
        return policy

    def volume(self, obj, path):
        fields = {'base_elevation', 'min_elevation', 'max_elevation', 'wall_min_height', 'roof_height'}
        self.object(obj, path, fields | {'roof', 'levels', 'height_type'}, ('roof',))
        values = {k: self.value(obj[k], path + (k,)) for k in fields if k in obj}
        roof = obj['roof']
        self.object(roof, path + ('roof',), {'type', 'reference_source', 'alignment', 'orientation', 'edge_index'}, ('type',))
        roof_values = {k: self.value(v, path + ('roof', k)) for k, v in roof.items()}
        wall = None
        if 'levels' in obj:
            levels = obj['levels']
            self.object(levels, path + ('levels',), {'wall'}, ('wall',))
            block = levels['wall']
            p = path + ('levels', 'wall')
            self.object(block, p, {'count', 'height', 'top_height', 'height_type', 'heights'}, ('count', 'height'))
            overrides = []
            self.array(block.get('heights', []), p + ('heights',))
            for i, override in enumerate(block.get('heights', [])):
                q = p + ('heights', i)
                self.object(override, q, {'index', 'height'}, ('index', 'height'))
                overrides.append((self.value(override['index'], q + ('index',)), self.value(override['height'], q + ('height',))))
            wall = LevelSpec(self.value(block['count'], p + ('count',)), self.value(block['height'], p + ('height',)),
                             self.value(block.get('top_height', 0), p + ('top_height',)), self.policy(block, p), tuple(overrides), p)
        if 'wall_min_height' not in values and wall is None:
            raise self.source.error(path, 'wall_min_height or a wall level model is required')
        return VolumeSpec(MappingProxyType(values), MappingProxyType(roof_values), wall, self.policy(obj, path))

    def check_conditions(self, program):
        def check(condition, chain=()):
            if condition is None:
                return
            if condition.kind == 'ref':
                name = condition.data
                if name not in program.conditions:
                    raise self.source.error(condition.path, f'Unknown condition {name!r}')
                if name in chain:
                    raise self.source.error(condition.path, 'Condition cycle: ' + ' -> '.join((*chain, name)), 'cyclic_reference')
                check(program.conditions[name], (*chain, name))
            elif condition.kind in ('and', 'or'):
                for child in condition.data:
                    check(child, chain)
            elif condition.kind == 'not':
                check(condition.data, chain)
        for name, condition in program.conditions.items():
            check(condition, (name,))
        for operation in program.walk_operations():
            check(operation.condition)
            if operation.kind == 'choose_one':
                for choice in operation.params:
                    check(choice.condition)
