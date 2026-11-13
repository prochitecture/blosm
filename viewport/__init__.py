"""Public native viewport API, version 1. Provider modules depend on this API."""
API_VERSION = 1


def open(context,surface):
    from .session import open as open_session
    return open_session(context,surface)


def current():
    from .session import active
    return active


def register():
    from .session import register as register_session
    register_session()


def unregister():
    from .session import unregister as unregister_session
    unregister_session()
