"""Credential-safe diagnostics shared by callers and the isolated helper."""
import re
import json
from urllib.parse import parse_qsl, quote, quote_plus, unquote, urlsplit, urlunsplit

_secret_name = re.compile(r'key|token|session|secret|password|signature|credential|authorization|cookie', re.I)
_url = re.compile(r'https?://[^\s<>"\']+', re.I)
_assignment = re.compile(
    r'''(?ix)(\b(?:api[_-]?key|key|access[_-]?token|token|session(?:[_-]?id)?|password|secret|signature|authorization)["']?\s*[:=]\s*)("[^"]*"|'[^']*'|[^\s&,}\]]+)''')


def redact_url(url):
    try:
        parts = urlsplit(url)
        # Never retain URL user-info, query values or fragments in diagnostics.
        host = parts.netloc.rsplit('@', 1)[-1]
        return urlunsplit((parts.scheme, host, parts.path,
                           '[redacted]' if parts.query else '', ''))
    except ValueError:
        return '[redacted URL]'


def redact_text(value, url='', headers=None):
    text = str(value)
    secrets = set()
    try:
        parts = urlsplit(url)
        secrets.update(v for k, v in parse_qsl(parts.query) if v and _secret_name.search(k))
        if parts.username:
            secrets.add(unquote(parts.username))
        if parts.password:
            secrets.add(unquote(parts.password))
    except ValueError:
        pass
    for key, value in (headers or {}).items():
        if value and _secret_name.search(key):
            secrets.add(value)
            if key.lower() == 'authorization':
                secrets.add(value.rsplit(' ', 1)[-1])
    text = _url.sub(lambda match: redact_url(match.group()), text)
    # An incomplete Authorization header can contribute an empty suffix. Never
    # replace the empty string: that inserts redaction markers between all chars.
    for secret in sorted(filter(None, secrets), key=len, reverse=True):
        for form in (secret, quote(secret, safe=''), quote_plus(secret), json.dumps(secret)[1:-1]):
            text = text.replace(form, '[redacted]')
    # This also covers echoed credentials not present in the request itself.
    return _assignment.sub(lambda match: match.group(1) + '"[redacted]"', text)


def redact_body(body, url, headers):
    def clean(item):
        if isinstance(item, dict):
            return {key: '[redacted]' if _secret_name.search(key) else clean(value)
                    for key, value in item.items()}
        if isinstance(item, list):
            return [clean(value) for value in item]
        return redact_text(item, url, headers) if isinstance(item, str) else item
    try:
        result = json.dumps(clean(json.loads(body)))
    except (ValueError, TypeError, RecursionError):
        result = redact_text(body, url, headers)
    # Redacting a short, repeated credential can expand a bounded body greatly.
    # Keep diagnostics within their IPC budget without truncating valid JSON.
    if len(result) > 65536:
        return json.dumps({'error': {'message': 'Provider error body omitted after redaction (size limit)'}})
    return result
