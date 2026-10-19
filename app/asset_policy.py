"""Storage policy shared by all Blender asset consumers; no Blender imports."""
from dataclasses import dataclass


@dataclass(frozen=True)
class AssetLoadPolicy:
    mode: str = 'APPEND'

    def __post_init__(self):
        if self.mode not in ('APPEND', 'AUTO', 'RETAINED'):
            raise ValueError('Unknown Blender asset storage policy: '+str(self.mode))

    def uses_link(self, spec, override=None):
        """Automatic linking requires published byte identity and certification.

        An explicit per-call choice wins. Mandatory retained mode still passes
        through the loader's certification, edition and durable-store checks;
        a missing certificate must never silently change it to embedding.
        """
        if override is not None:
            if not isinstance(override, bool):
                raise TypeError('Asset link override must be True, False or None')
            return override
        return self.mode == 'RETAINED' or (
            self.mode == 'AUTO' and bool(spec.sha256) and spec.self_contained)
