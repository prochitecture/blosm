"""Preference-first token selection with one lazy file read per Blender session.

Only process memory holds the fallback. Do not log it or save it to preferences
or scene data. Normal add-on disable/re-enable retains this module in sys.modules
and therefore retains the cache.
"""
from pathlib import Path

TOKEN_FILE = Path(__file__).with_name('mapillary_token.txt')
_file_loaded = False
_file_token = ''
_file_error = None


def get_token(preference_token=''):
    """Check preferences each call; access the fallback file at most once.

    Missing/empty/unreadable outcomes are cached too. Preferences can override
    the cache immediately; restart Blender to pick up later file changes.
    """
    global _file_loaded, _file_token, _file_error
    preferred = preference_token.strip()
    if preferred:
        return preferred
    if not _file_loaded:
        _file_loaded = True
        try:
            _file_token = TOKEN_FILE.read_text(encoding='utf-8-sig').strip()
        except FileNotFoundError:
            pass
        except (OSError, UnicodeError):
            # Decode exceptions may contain secret bytes. Cache only a constant
            # diagnostic and avoid retrying on every import.
            _file_error = ('Unable to read mapillary/mapillary_token.txt. '
                           'Set a Mapillary client access token in Blosm preferences.')
    if _file_error:
        raise RuntimeError(_file_error) from None
    return _file_token
