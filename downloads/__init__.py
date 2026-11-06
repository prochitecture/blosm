"""Reusable in-memory HTTP downloads; no Blender dependencies."""
from dataclasses import dataclass, field


@dataclass(frozen=True)
class Request:
    url: str
    headers: dict = field(default_factory=dict)
    decode_gzip: bool = False
    verify_tls: bool = True
