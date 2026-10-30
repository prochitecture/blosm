"""Versioned deterministic sampling using arithmetic reproducible in Geometry Nodes.

All mixer intermediates are exact integers below 2**24, including in float32.
This is a procedural sampling policy, not a cryptographic generator. Stable
source identities and authored occurrence paths, never mesh indices, name draws.
"""
import hashlib
import json

POLICY = 'padl-mix20-v1'
MODULUS = 1048576
ROUNDS = (123457, 456791, 789221, 345679)


def stable_key(value):
    data = json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(',', ':')).encode('utf-8')
    return int.from_bytes(hashlib.sha256(data).digest()[:4], 'little') % MODULUS


def unit_from_keys(seed, owner_key, draw_key):
    state = (seed % MODULUS + (seed // MODULUS) * 1013 + owner_key + draw_key) % MODULUS
    for offset in ROUNDS:
        state = ((state % 1024) * 7919 + (state // 1024) * 1049 + offset) % MODULUS
    return (state + 0.5) / MODULUS


def unit(seed, scope, owner, occurrence):
    return unit_from_keys(seed, stable_key((scope, owner)), stable_key(occurrence))
