"""Capture-date and photo-model predicates, independent of Blender and HTTP."""
from dataclasses import dataclass
from datetime import date
from decimal import Decimal, InvalidOperation
import re

DAY_MS = 86400000
EPOCH = date(1970, 1, 1)
MIN_MS = (date.min - EPOCH).days * DAY_MS
MAX_MS = ((date.max - EPOCH).days + 1) * DAY_MS
CATEGORIES = ('ordinary', 'panoramic', 'fisheye', 'other')


def date_bound(text, label):
    text = text.strip()
    if not text:
        return None
    if not re.fullmatch(r'[0-9]{4}-[0-9]{2}-[0-9]{2}', text):
        raise ValueError('%s date must be YYYY-MM-DD, or blank' % label)
    try:
        value = date.fromisoformat(text)
    except ValueError:
        raise ValueError('%s date is not a valid calendar date' % label) from None
    return (value - EPOCH).days * DAY_MS


def capture_timestamp(value):
    """Read milliseconds, including older string/float custom properties.

    No local timezone conversion: the date controls denote UTC calendar days.
    Missing, malformed and out-of-calendar-range timestamps are undated.
    """
    if isinstance(value, bool) or not isinstance(value, (str, int, float)):
        return None
    try:
        stamp = Decimal(str(value))
        if stamp.is_finite() and MIN_MS <= stamp < MAX_MS:
            return stamp
    except InvalidOperation:
        pass
    return None


def photo_category(model):
    model = model.strip().lower() if isinstance(model, str) else ''
    return {'perspective': 'ordinary', 'spherical': 'panoramic',
            'equirectangular': 'panoramic', 'fisheye': 'fisheye'}.get(model, 'other')


@dataclass(frozen=True)
class CameraFilter:
    start: int | None
    end: int | None
    categories: frozenset
    include_undated: bool

    @classmethod
    def build(cls, date_from='', date_to='', categories=CATEGORIES, include_undated=False):
        start, last = date_bound(date_from, 'From'), date_bound(date_to, 'To')
        if start is not None and last is not None and start > last:
            raise ValueError('From date must be on or before To date')
        # Exclusive upper bound includes the entire To day, even 9999-12-31.
        return cls(start, None if last is None else last + DAY_MS,
                   frozenset(categories), bool(include_undated))

    def matches(self, captured_at, model):
        if photo_category(model) not in self.categories:
            return False
        if self.start is None and self.end is None:
            return True  # Undated cameras are included without a date range.
        stamp = capture_timestamp(captured_at)
        if stamp is None:
            return self.include_undated
        return ((self.start is None or stamp >= self.start) and
                (self.end is None or stamp < self.end))
