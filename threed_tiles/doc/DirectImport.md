# Direct glTF import

The 3D Tiles renderer calls Blender's glTF implementation through
threed_tiles/gltf_import.py. Downloaded GLB and B3DM tiles are imported from
memory. The renderer does not write a temporary tile or invoke
bpy.ops.import_scene.gltf.

The implementation is also used by the optional progressive importer described
in [ProgressiveImport.md](../../progressive/ProgressiveImport.md). Shared
[concurrent downloads](../../downloads/README.md) work with both import modes.

## Entry point and version compatibility

import_gltf(filepath, *, content=None, report=None) accepts either a file path
or a bytes-like object. When content is provided, the path supplies the source
identity and the directory for external buffers/images, without requiring that
the tile itself exist. Self-contained GLB content requires no tile file.

The function returns attribution and importer messages. It does not retain the
importer's decoded buffers. Optional report accepts an operator-style
report({level}, message) callback; exceptions propagate to the caller.

The renderer receives the outer Blosm operator's report callback and passes it
to the adapter for GLB and B3DM imports. Raised errors propagate to the manager
for per-tile error accounting without individual UI reports. Blosm continues
with subsequent tiles, prints each failure with its URI in the System Console,
and reports the aggregate result once in the UI. Unsupported extensions across
many tiles therefore produce one final error summary, rather than repeated
popups. A partial import reports its successful tile count and error count.

Buffered glTF warnings and nonfatal errors are forwarded with their original
severity only after a successful tile import, matching the glTF operator's
success path. Failed-import diagnostics remain in the console. Logger cleanup
runs on both success and failure. Python calls through bpy.ops.blosm.import_data()
still raise RuntimeError when the outer operator reports an ERROR, such as when
no tiles were imported.

Version dispatch follows the existing coordinate-conversion patching:

| Blender version | Scene creation implementation |
| --- | --- |
| 3.6 and 4.0–4.2 | io_scene_gltf2.blender.imp.gltf2_blender_gltf.BlenderGlTF.create |
| 4.3 onward | io_scene_gltf2.blender.imp.blender_gltf.BlenderGlTF.create |

The installed operator's RNA defaults supply settings without invoking the
operator. Enabled glTF import extensions are discovered as before. Old and new
logger cleanup APIs are handled separately. Release's declared minimum Blender
version remains 4.5; the legacy adapter was also checked on Blender 3.6.4.

The adapter explicitly updates the view layer after scene creation. Without
this update, immediate reads of new objects' world matrices can lose the
geographic offset or RTC placement. No pre-import update is needed; import
parity was checked without it on Blender 4.5, 5.0, 5.1 and 5.2.

## Renderer behavior

- GLB downloads go directly to the adapter as bytes.
- B3DM is decoded with the existing bundled parser. Its GLB array goes directly
  to the adapter, while RTC_CENTER is applied using the existing renderer logic.
  The temporary RTC conversion switch is restored on exceptions.
- Attribution is read from parsed glTF asset metadata, preserving spaces,
  punctuation and Unicode. Existing semicolon-separated aggregation remains.
- Root and nested JSON manifests are downloaded and parsed in memory on each
  import. GLB and B3DM tiles are also downloaded on each import. There are no
  dataset cache settings, cache reads or cache writes. Existing cache files are
  ignored; 3D Tiles import does not require or create a data directory.
- The renderer supplies a tile filename only for glTF source naming. The
  adapter's file-input mode remains available for comparison benchmarks;
  production tile imports always supply downloaded content.
- The outer Blosm import operator still declares REGISTER and UNDO. Individual
  glTF operations no longer create an operator boundary. Like Blender's native
  importer, the adapter does not roll back objects created before an import fails.

Scene creation still runs synchronously on Blender's main thread. This change
does not make a single glTF import interruptible or safe in a worker thread.
The blocking path owns coordinate patches for the whole import. The progressive
adapter applies and restores them within each synchronous tile step, before
returning control to Blender.

## Validation and benchmark

Run these scripts in factory-startup background Blender, with all generated
outputs below an external artifact directory:

~~~text
blender --background --factory-startup --disable-autoexec --python-exit-code 1 \
  --python tests/check_direct_gltf_import.py -- \
  --artifacts <artifact-directory>/checks --helsinki <directory-of-b3dm-fixtures>

blender --background --factory-startup --disable-autoexec --python-exit-code 1 \
  --python tests/check_direct_gltf_operator.py -- \
  --artifacts <artifact-directory>/operator-check

blender --background --factory-startup --disable-autoexec --python-exit-code 1 \
  --python tests/benchmark_direct_gltf_import.py -- \
  --artifacts <artifact-directory>/benchmark \
  --helsinki <directory-of-b3dm-fixtures> --repeats 20
~~~

These are argument layouts; shell continuation syntax differs on Windows.
Do not run benchmark processes concurrently. The baseline renderer is loaded
from fixed historical commit 4adff00270f4492667f3bf1305d0b338f3cd1268,
before direct importing was introduced. Its revision is recorded in the results.

Correctness checks cover operator/file/byte equivalence, topology, UVs, vertex
colors, packed images, matrices, attribution, external resources, enabled user
extensions, invalid input, logger cleanup, RTC with/without/zero center, joining,
and patch restoration after failure. A separate check exercises the registered
Blosm operator with root and nested JSON manifests and mixed GLB/B3DM tiles,
including repeated downloads and continuation after a failed tile. Guards reject
local dataset reads/writes and any invocation of the glTF import operator in the
production paths. The operator also works with a nonexistent data directory.

The checks passed on Blender 4.5.0, 5.0.0, 5.1.0 and 5.2.0 LTS, including three
real Helsinki tiles. Adapter-only checks also passed on Blender 3.6.4.

The benchmark uses three warm-ups per path, twenty measured rounds, and a
deterministic shuffled path order. It compares:

1. The glTF operator with an existing GLB.
2. The direct adapter with that same existing GLB.
3. Bytes written to a temporary GLB, imported with the operator, then deleted.
4. The same bytes imported directly.
5. The complete old and new B3DM renderer paths, including conversion, RTC and
   finalization.

Network, fixture loading and scene cleanup are excluded. Each trial starts with
an empty scene. Extracted GLB samples come from the same B3DM data, so they are
not a separate dataset. Timing uses perf_counter wall time. The JSON output
includes every measurement, medians, interquartile ranges and source hashes.

On this Windows machine, the initial full B3DM renderer measurements showed
8.7–31.7% lower local import time on Blender 5.2 and 14.9–40.1% on Blender 4.5
across the three samples. The reduction comes primarily from avoiding the
temporary-file round trip: direct imports of existing files alone did not show
a consistent performance improvement. These results do not predict network
speed, viewport responsiveness, or performance in an already large scene.

Raw fixtures, logs, verification reports and the detailed benchmark report are
under C:/dev/workspace/artifacts/direct_gltf_import, outside Git.

## Dataset attribution

Source: [City of Helsinki, Helsinki 2024 3D mesh](https://kartta.hel.fi/3d/mesh/Helsinki_2024/tileset.json).
The user supplied this dataset under CC BY-SA 4.0. The benchmark uses only
1/mesh_2.b3dm, 1/mesh_27.b3dm, and 1/mesh_5.b3dm, plus GLBs extracted from those
files for comparison. Fixture files and derivatives remain in artifacts;
they are not distributed with the add-on.
