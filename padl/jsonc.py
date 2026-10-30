"""Strict JSON with line/block comments, duplicate detection and JSON-path locations."""
import bisect
import json
import math

from .errors import SourceMap


def strip_comments(text, filename='<string>'):
    chars = list(text)
    i = 0
    quoted = False
    while i < len(text):
        if quoted:
            if text[i] == '\\':
                i += 2
                continue
            if text[i] == '"':
                quoted = False
        elif text[i] == '"':
            quoted = True
        elif text.startswith('//', i):
            end = text.find('\n', i)
            end = len(text) if end < 0 else end
            chars[i:end] = [' '] * (end - i)
            i = end
            continue
        elif text.startswith('/*', i):
            end = text.find('*/', i + 2)
            if end < 0:
                source = SourceMap(filename, {(): (text.count('\n', 0, i) + 1,
                                                   i - text.rfind('\n', 0, i))})
                raise source.error((), 'Unterminated block comment', 'json_syntax')
            end += 2
            chars[i:end] = [c if c in '\r\n' else ' ' for c in text[i:end]]
            i = end
            continue
        i += 1
    return ''.join(chars)


class JsonReader:
    def __init__(self, text, filename='<string>'):
        self.text = strip_comments(text.lstrip('\ufeff'), filename)
        self.source = SourceMap(filename)
        self.newlines = [-1] + [i for i, c in enumerate(self.text) if c == '\n']
        self.decoder = json.JSONDecoder()
        self.i = 0

    def whitespace(self):
        while self.i < len(self.text) and self.text[self.i] in ' \t\r\n':
            self.i += 1

    def mark(self, path):
        line = bisect.bisect_left(self.newlines, self.i)
        self.source.locations[path] = (line, self.i - self.newlines[line - 1])

    def fail(self, path, message):
        self.mark(path)
        raise self.source.error(path, message, 'json_syntax')

    def read(self):
        result = self.value(())
        self.whitespace()
        if self.i != len(self.text):
            self.fail((), 'Unexpected content after the JSON document')
        return result, self.source

    def value(self, path):
        if len(path) > 128:
            self.fail(path, 'Maximum document nesting depth (128) exceeded')
        self.whitespace()
        self.mark(path)
        c = self.text[self.i:self.i + 1]
        if c in ('{', '['):
            is_object = c == '{'
            closing = '}' if is_object else ']'
            result = {} if is_object else []
            self.i += 1
            self.whitespace()
            if self.text[self.i:self.i + 1] == closing:
                self.i += 1
                return result
            while True:
                self.whitespace()
                if is_object:
                    if self.text[self.i:self.i + 1] != '"':
                        self.fail(path, 'Expected an object key')
                    key = self.scalar(path)
                    if key in result:
                        self.fail(path + (key,), f'Duplicate object key {key!r}')
                    self.whitespace()
                    if self.text[self.i:self.i + 1] != ':':
                        self.fail(path + (key,), 'Expected a colon')
                    self.i += 1
                    result[key] = self.value(path + (key,))
                else:
                    result.append(self.value(path + (len(result),)))
                self.whitespace()
                c = self.text[self.i:self.i + 1]
                self.i += 1
                if c == closing:
                    return result
                if c != ',':
                    self.fail(path, f'Expected a comma or {closing}')
        return self.scalar(path)

    def scalar(self, path):
        try:
            value, self.i = self.decoder.raw_decode(self.text, self.i)
        except json.JSONDecodeError as exc:
            self.i = exc.pos
            self.fail(path, exc.msg)
        if isinstance(value, float) and not math.isfinite(value):
            self.fail(path, 'JSON numbers must be finite')
        return value
