# Conforming prescribed separator boundaries

`conforming_boundary_refinement.refine` inserts an independently qualified
grass disk's additional boundary samples into the existing pavement mesh.
Inputs are the mesh, authoritative boundary references and group, pattern
library, and a dictionary from hole index to counterclockwise XY rim.
It returns the mesh, new references, original vertex map and construction
metadata. Every original point must occur byte-exact in source order, and
every new point must lie on its original finite curb chord. Outer rims and
incoming profile chords cannot be split.

Odd face parity propagates across internal dual edges. Positive cached quad
templates cover adjacent even side counts as well as matching opposite parity.
Unequal fixed corner edges can make a transition quad impossible to repair;
bounded local collars separate those edges before mapping a direct-quad core.
The ordinary shape correction, edge rotations and anisotropic refinement
qualify each such original face with its complete rim fixed. Construction uses
a local numerical frame and restores all authoritative points byte-exact.
Heights and face-corner UV seams retain their original bilinear source map.

The helper preserves the complete finite polygon and returns positive,
conforming quads. Stored-coordinate shape, height, slope, drainage and joins
still require the ordinary final checks. Metadata identifies the touched
original faces and one-ring original vertex IDs, allowing callers to retain
far fixed guides during subsequent correction. Cell, insertion and per-face
atlas budgets are finite; failure does not mutate the input mesh.

Captured original Bundle 903 qualification inserted 192 boundary points and
172 internal parity edges. The resulting 15,049 pavement quads had maximum
double-coordinate edge ratio 7.8943; actual float32 coordinates retained edge
ratio 7.9666 and minimum corner sine 0.0388. All original pavement points and
both requested hole rims were byte-exact. Actual stored-coordinate checks run
for both fast templates and local graded faces, with a tighter construction
target used only where storage needs numerical room. The final limits remain
edge ratio 8 and corner sine 0.03. Development captures and subsequent surface
reports belong under the street task artifacts directory.
