"""Mapillary spherical projection (SfM +X right, +Y down, +Z forward).

Transform.ts rotates the *2D image* for EXIF, after spherical projection.
Undo the EXIF rotation in the imported camera's 3D basis before projection;
then rotate UVs once. Rotating a spherical bearing alone is not equivalent.
"""
import math


def raw_from_local(direction,orientation=1):
    x,y,z=direction
    # Blender camera -> basic SfM, then inverse of Transform's Rz(EXIF).
    x,y,z=x,-y,-z
    if orientation==3: return -x,-y,z
    if orientation==6: return y,-x,z
    if orientation==8: return -y,x,z
    if orientation!=1: raise ValueError('Mirrored EXIF panoramas are not supported')
    return x,y,z


def image_uv(direction,orientation=1):
    """Return UV in the decoded upright image, with Blender bottom-left origin."""
    x,y,z=raw_from_local(direction,orientation)
    u=(.5+math.atan2(x,z)/(2*math.pi))%1.
    v=.5+math.atan2(y,math.hypot(x,z))/math.pi
    if orientation==3: u,v=1-u,1-v
    elif orientation==6: u,v=1-v,u
    elif orientation==8: u,v=v,1-u
    return u,1-v


def valid_image_size(width,height,orientation):
    expected=.5 if orientation in (6,8) else 2.
    return width>0 and height>0 and abs(width/height-expected)<.02
