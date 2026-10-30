"""Diagnostics shared by the PADL frontend and geometry backends."""
from dataclasses import dataclass


def pointer(path):
    return '/' + '/'.join(str(part).replace('~', '~0').replace('/', '~1') for part in path)


@dataclass(frozen=True)
class SourceLocation:
    filename: str
    path: tuple = ()
    line: int = 1
    column: int = 1

    def __str__(self):
        return f'{self.filename}:{self.line}:{self.column} {pointer(self.path)}'


class PadlError(ValueError):
    def __init__(self, message, *, code='invalid_document', location=None):
        self.message = message
        self.code = code
        self.location = location
        super().__init__(f'{location}: {message}' if location else message)


class SourceMap:
    def __init__(self, filename, locations=None):
        self.filename = str(filename)
        self.locations = locations or {}

    def at(self, path=()):
        parent = tuple(path)
        while parent not in self.locations and parent:
            parent = parent[:-1]
        line, column = self.locations.get(parent, (1, 1))
        return SourceLocation(self.filename, tuple(path), line, column)

    def error(self, path, message, code='invalid_document'):
        return PadlError(message, code=code, location=self.at(path))
