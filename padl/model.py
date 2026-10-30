"""Backend-independent operation model. Importing this module does not load Blender."""
from dataclasses import dataclass, field
from typing import Any, Mapping

from .errors import SourceMap


@dataclass(frozen=True)
class Value:
    kind: str
    data: Any
    path: tuple = ()


@dataclass(frozen=True)
class Condition:
    kind: str
    data: Any
    path: tuple = ()


@dataclass(frozen=True)
class LevelSpec:
    count: Value
    height: Value
    top_height: Value
    height_type: str = 'fixed'
    overrides: tuple = ()
    path: tuple = ()


@dataclass(frozen=True)
class VolumeSpec:
    values: Mapping[str, Value]
    roof: Mapping[str, Value]
    wall: LevelSpec | None = None
    height_type: str = 'fixed'


@dataclass(frozen=True)
class Choice:
    name: str
    condition: Condition | None
    operations: tuple
    path: tuple


@dataclass(frozen=True)
class Operation:
    kind: str
    name: str
    params: Any
    condition: Condition | None
    path: tuple


@dataclass(frozen=True)
class Program:
    parameters: Mapping[str, Value]
    conditions: Mapping[str, Condition]
    operations: tuple[Operation, ...]
    source: SourceMap
    fingerprint: str
    schema_version: str
    context: Mapping = field(default_factory=dict)
    document: Mapping = field(default_factory=dict, repr=False)

    def walk_operations(self):
        def walk(items):
            for op in items:
                yield op
                if op.kind == 'choose_one':
                    for choice in op.params:
                        yield from walk(choice.operations)
        return walk(self.operations)
