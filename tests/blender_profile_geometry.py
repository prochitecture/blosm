"""Legacy roof adapters continue to use the shared profile kernel."""
from pathlib import Path
import sys,math
from types import SimpleNamespace
import bpy,addon_utils
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
addon_utils.enable('blosm',default_set=True)
from blosm.building.roof.profile import RoofProfile,ProfiledVert,Slot,gabledRoof,roundRoof,gambrelRoof,saltboxRoof
from blosm.building.roof import PolygonOLD
from blosm.padl.profile_geometry import ProfiledRoofGeometry,ProfiledVert as SharedVertex,Slot as SharedSlot
from blosm.padl.geometry import signed_area
assert issubclass(RoofProfile,ProfiledRoofGeometry)
assert ProfiledVert is SharedVertex and Slot is SharedSlot
PolygonOLD.straightAngleTan = math.tan(math.radians(abs(180-bpy.context.scene.blosm.straightAngleThreshold)))
rings=[((0,0),(12,0),(12,8),(0,8)),((0,0),(12,0),(12,4),(7,4),(7,8),(0,8))]
for profile in (gabledRoof,roundRoof,gambrelRoof,saltboxRoof):
    roof=RoofProfile(profile)
    for ring in rings:
        for orientation in ('along','across'):
            roof.init(SimpleNamespace(tags={'height':'30','roof:height':'5','min_height':'2','roof:orientation':orientation}),ring,None,None)
            assert roof.valid and roof.make(None)
            assert roof.roofIndices and roof.wallIndices
            area=sum(abs(signed_area([roof.verts[i] for i in face])) for face in roof.roofIndices)
            assert abs(area-abs(signed_area(ring))) < .001,(area,ring)
            assert min(v.z for v in roof.verts)==2
            assert abs(max(v.z for v in roof.verts)-30)<.001
print('LEGACY_PROFILE_ADAPTERS_OK 16',flush=True)
