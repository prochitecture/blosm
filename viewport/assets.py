"""Reference ownership, independent of a particular provider or HTTP service."""
from dataclasses import dataclass


@dataclass(frozen=True)
class AssetTicket:
    session_id: int
    generation: int
    key: str


class ImageRef:
    """Borrow a Blender image. Releasing GPU references never deletes its ID."""
    def __init__(self,image):
        import gpu
        self.image = image
        self.size = tuple(image.size)
        self.texture = gpu.texture.from_image(image)
        self.premultiplied = image.alpha_mode == 'PREMUL'

    def release(self):
        self.texture = self.image = None
