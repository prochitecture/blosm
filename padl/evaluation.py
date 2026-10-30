"""Lazy scoped evaluation and extensible operation execution, without Blender."""
import copy
import math
from dataclasses import dataclass, field
from typing import Mapping

from . import randomness
from .model import Value


def is_number(value):
    try:
        return type(value) in (int, float) and math.isfinite(value)
    except OverflowError:
        return False


@dataclass(frozen=True)
class FootprintContext:
    identity: str
    attributes: Mapping = field(default_factory=dict)
    edge_identity: str | None = None
    edge_attributes: Mapping = field(default_factory=dict)


@dataclass(frozen=True)
class EvaluationResult:
    volumes: tuple
    parameters: Mapping
    samples: Mapping


class ValueResolver:
    def __init__(self, program, *, seed=0, global_inputs=None):
        if type(seed) is not int or not -(2**31) <= seed < 2**31:
            raise program.source.error(('context', 'seed'), 'Seed must be a signed 32-bit integer', 'invalid_input')
        self.program = program
        self.seed = seed
        self.global_inputs = copy.deepcopy(dict(global_inputs or {}))
        self.cache = {}
        self.samples = {}
        self.stack = []

    def fail(self, value, message, code='invalid_value'):
        raise self.program.source.error(value.path, message, code)

    def numeric(self, result, value):
        if not is_number(result):
            self.fail(value, f'Expected a finite number, got {result!r}')
        return result

    def owner(self, scope, context, value):
        if scope == 'global':
            return 'global', self.global_inputs
        if not context or not isinstance(context.identity, str) or not context.identity:
            self.fail(value, f'{scope} requires a stable footprint identity', 'missing_context')
        if scope == 'footprint':
            return context.identity, context.attributes
        if context.edge_identity is None:
            self.fail(value, 'footprint_edge requires a stable originating-edge identity', 'missing_context')
        return (context.identity, context.edge_identity), context.edge_attributes

    def parameter(self, value, context):
        scope, name = value.data
        ranks = {'global': 0, 'footprint': 1, 'footprint_edge': 2}
        if self.stack and ranks[scope] > ranks[self.stack[-1][0]]:
            self.fail(value, f'{self.stack[-1][0]} parameter {self.stack[-1][2]!r} cannot depend on {scope} parameter {name!r}', 'invalid_scope_dependency')
        owner, inputs = self.owner(scope, context, value)
        key = (scope, owner, name)
        if key in self.cache:
            return self.cache[key]
        if name in inputs:
            result = inputs[name]
            if type(result) not in (str, bool, int, float) or (type(result) in (int, float) and not is_number(result)):
                self.fail(value, f'Supplied parameter {name!r} must be a finite scalar', 'invalid_input')
        else:
            if name not in self.program.parameters:
                self.fail(value, f'Unresolved {scope} parameter {name!r}', 'unresolved_parameter')
            if key in self.stack:
                chain = ' -> '.join(f'{s}:{o}:{n}' for s, o, n in (*self.stack, key))
                self.fail(value, f'Parameter cycle: {chain}', 'cyclic_reference')
            self.stack.append(key)
            try:
                result = self.resolve(self.program.parameters[name], context)
            finally:
                self.stack.pop()
        self.cache[key] = result
        return result

    def draw(self, value, context):
        if self.stack:
            scope, owner, name = self.stack[-1]
            occurrence = ('parameter', name, value.path)
        else:
            scope = 'footprint' if context is not None else 'global'
            owner, _ = self.owner(scope, context, value)
            occurrence = ('inline', value.path)
        key = (scope, owner, occurrence)
        if key not in self.samples:
            self.samples[key] = randomness.unit(self.seed, scope, owner, occurrence)
        return self.samples[key]

    def resolve(self, value, context=None):
        kind, data = value.kind, value.data
        if kind == 'literal':
            return data
        if kind == 'parameter':
            return self.parameter(value, context)
        if kind == 'expression':
            op, children = data
            args = [self.numeric(self.resolve(child, context), child) for child in children]
            if op == 'divide' and args[1] == 0:
                self.fail(value, 'Division by zero', 'division_by_zero')
            funcs = {'add': lambda a, b: a + b, 'subtract': lambda a, b: a - b,
                     'multiply': lambda a, b: a * b, 'divide': lambda a, b: a / b,
                     'negate': lambda a: -a, 'min': min, 'max': max,
                     'abs': abs, 'floor': math.floor, 'ceil': math.ceil}
            try:
                result = funcs[op](*args)
            except (OverflowError, ValueError) as exc:
                self.fail(value, str(exc), 'numeric_overflow')
            return self.numeric(result, value)
        if kind == 'uniform':
            low, high = (self.numeric(self.resolve(v, context), v) for v in data)
            if high < low:
                self.fail(value, 'Uniform maximum must not be below its minimum')
            return self.numeric(low + (high - low) * self.draw(value, context), value)
        if kind == 'weighted_choice':
            weights = [self.numeric(self.resolve(w, context), w) for _, w in data]
            total = sum(weights)
            if any(w < 0 for w in weights) or not is_number(total) or total <= 0:
                self.fail(value, 'Choice weights must be nonnegative with a finite positive sum')
            threshold = self.draw(value, context) * total
            cumulative = 0
            for (choice, _), weight in zip(data, weights):
                cumulative += weight
                if threshold < cumulative:
                    return self.resolve(choice, context)
            self.fail(value, 'Weighted choice did not select a value')
        self.fail(value, f'Unsupported value {kind}', 'unsupported_feature')

    def condition(self, condition, context=None):
        if condition is None:
            return True
        kind, data = condition.kind, condition.data
        if kind == 'ref':
            return self.condition(self.program.conditions[data], context)
        if kind == 'and':
            return all(self.condition(c, context) for c in data)
        if kind == 'or':
            return any(self.condition(c, context) for c in data)
        if kind == 'not':
            return not self.condition(data, context)
        op, children = data
        args = [self.resolve(v, context) for v in children]
        if op in ('eq', 'ne'):
            a, b = args
            if type(a) != type(b) and not (is_number(a) and is_number(b)):
                raise self.program.source.error(condition.path, f'Cannot compare {type(a).__name__} and {type(b).__name__}', 'invalid_value')
            return (a == b) if op == 'eq' else (a != b)
        for result, child in zip(args, children):
            self.numeric(result, child)
        if op == 'between':
            if args[1] > args[2]:
                raise self.program.source.error(condition.path, 'Comparison minimum exceeds maximum')
            return args[1] <= args[0] <= args[2]
        return {'lt': lambda a, b: a < b, 'lte': lambda a, b: a <= b,
                'gt': lambda a, b: a > b, 'gte': lambda a, b: a >= b}[op](*args)


class Evaluator:
    """One document execution; each evaluate() starts fresh input-dependent state.

    Stateless, identity-derived draws preserve unchanged random preferences
    across edits. Allocated heights never replace resolved declared values.
    Operation handlers receive ``(evaluator, operation, resolver, context)``.
    """
    def __init__(self, program, *, seed=None, global_inputs=None, handlers=None):
        self.program = program
        self.seed = program.context.get('seed', 0) if seed is None else seed
        self.global_inputs = copy.deepcopy(global_inputs or {})
        self.handlers = {'choose_one': Evaluator.choose_one, 'generate_volume': Evaluator.volume}
        self.handlers.update(handlers or {})

    def evaluate(self, context):
        context = FootprintContext(context.identity, copy.deepcopy(context.attributes),
                                   context.edge_identity, copy.deepcopy(context.edge_attributes))
        resolver = ValueResolver(self.program, seed=self.seed, global_inputs=self.global_inputs)
        volumes = tuple(self.operations(self.program.operations, resolver, context))
        return EvaluationResult(volumes, dict(resolver.cache), dict(resolver.samples))

    def operations(self, operations, resolver, context):
        for operation in operations:
            if not resolver.condition(operation.condition, context):
                continue
            handler = self.handlers.get(operation.kind)
            if handler is None:
                raise self.program.source.error(operation.path, f'No Python handler for {operation.kind}', 'unsupported_operation')
            yield from handler(self, operation, resolver, context)

    def choose_one(self, operation, resolver, context):
        for choice in operation.params:
            if resolver.condition(choice.condition, context):
                yield from self.operations(choice.operations, resolver, context)
                break

    def volume(self, operation, resolver, context):
        from .heights import HeightResolver
        yield HeightResolver(self.program.source).resolve(operation, resolver, context)
