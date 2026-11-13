"""Compatible optional producers for way preparation and rendering."""


def producer_for(target, method):
    """Keep newer synchronous overrides ahead of inherited step producers.

    An extension that overrides only the blocking method retains that behavior
    as one atomic operation. It can opt into cooperative work by providing its
    own corresponding ``*_steps`` method, matching the shared OSM hook policy.
    """
    hook = method + '_steps'
    producer = getattr(target, hook, None)
    if producer is None:
        return None
    mro = type(target).__mro__
    method_owner = next((base for base in mro if method in base.__dict__), None)
    hook_owner = next((base for base in mro if hook in base.__dict__), None)
    if method_owner and hook_owner and mro.index(method_owner) < mro.index(hook_owner):
        return None
    return producer
