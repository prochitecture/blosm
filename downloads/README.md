# Shared concurrent downloads

`downloads` contains no Blender imports. `DownloadTransport` runs one ordinary
Python helper with a bounded `ThreadPoolExecutor`; all Blender work stays on the
main thread. The helper is reused throughout an import and accepts requests as
new resources are discovered. Neither requests nor results create dataset files.

## Interface

- `Request(url, headers={}, decode_gzip=False, verify_tls=True)` describes an HTTP
  request. Headers and options are copied into the wire message on submission.
  Raw responses are the default so an intrinsically compressed resource, such as
  an HGT gzip file, can be passed unchanged to its consumer.
- `DownloadTransport(max_requests=4, max_bytes=256 MiB, max_buffer_bytes=None)`
  starts the helper. `submit(request)` returns a monotonically increasing identity.
  `poll()` returns `None` or `(identity, content, exception)` in completion order.
  A reply is delivered once. Requests never share an output file or Blender state.
- `close()` closes sockets, discards buffered replies and terminates the entire
  helper, including its worker threads. `reap()` performs nonblocking cleanup.
  No helper thread runs inside Blender, including in blocking import mode.

The maximum is eight outstanding requests per transport. Blosm currently allows
only one active import, so this is also its active import's total limit, shared
by manifest and content requests. Independent future services must coordinate
limits if they allow multiple simultaneous operations; this is not a global
provider quota service.

## Scheduling and memory

`progressive.core.Job` can drive serial requests and concurrent consumers:

```python
from functools import partial
from .downloads import Request
from .progressive.core import Download, Drain

def steps():
    manifest = yield Request(root_url)
    for item in discover(manifest):
        yield Download(Request(item.url), partial(import_item, item))
    yield Drain()
    finish_import()
```

`import_item(item, content, exception)` runs on the caller's main thread, never
in a worker. A normal `Request` suspends discovery until its result is returned
or its exception is thrown into the generator. `Download` permits further
submissions while earlier downloads run. `Drain` waits for every deferred
consumer before finalization. `None` remains a cooperative traversal checkpoint.

Each scheduler step performs bounded IPC/traversal work and at most one deferred
consumer. Paused polling receives results but performs no scene changes or new
submissions. Active downloads **and unconsumed results** occupy the same bounded
window; discovery stops when it fills. Reserving `max_bytes` per request enforces
a conservative payload budget even when content length is unknown or compressed.
An explicit smaller `max_buffer_bytes` reduces effective concurrency. The default
four slots reserve up to 1 GiB of response payloads (256 MiB each). This is not a
process-RSS guarantee: socket buffers, serialization copies, compressed data and
the currently decoded/imported Blender tile add overhead. There is no unlimited
queue of downloaded tiles. Consumers must release payloads after use.

The same Job can be stepped from Blender's timer or synchronously by a blocking
runner. Both use the identical transport and generator; only the progressive
controller adds UI yielding, context preservation and grouped undo checkpoints.

## 3D Tiles integration

Concurrent downloads are enabled by default, independently of progressive import.
The panel exposes a toggle and a simultaneous-download limit (default four).
Disabling the toggle retains the previous sequential path. Setting the limit to
one exercises the concurrent implementation with one outstanding request.

Manifest traversal remains depth-first and manifest dependencies are awaited in
order. Selected GLB/B3DM requests are deferred; complete URLs, paths and accumulated
transforms are captured before submission. Thus changes to session query state in
a later manifest cannot affect an earlier request. Already discovered content can
continue downloading while a child manifest is fetched. Independent manifests are
not fetched in parallel in this first implementation, preserving existing session
semantics. Results may import out of order, and future downloads overlap the
main-thread glTF import. Each tile's attribution/error accounting is retained.

Finalization waits for discovery and all submitted content, even if traversal
encounters an error after submitting some tiles. Cancellation discards outstanding
work and retains already committed progressive tiles. Final joining, result
popups, destination validation and grouped undo retain the progressive controller's
existing behavior. GLB/B3DM bytes go directly to glTF; caching remains excluded.

## Failures and lifecycle

The helper retries transient connection failures and HTTP 429/500/502/503/504 up
to twice, with exponential backoff/jitter and Retry-After support. Retry-After
values above 30 seconds are returned as errors rather than retried prematurely.
Permanent HTTP failures and invalid/truncated/oversized payloads are errors,
never silently interpreted as empty tiles. The socket operation timeout is 30
seconds; the caller enforces a 90-second request deadline and 15-second startup
deadline. A failed helper/protocol/deadline stops the whole job and closes it.

The shared transport verifies TLS by default. The 3D Tiles adapter explicitly
retains its pre-existing certificate-verification bypass for compatibility;
changing that existing provider policy is separate work. Helper diagnostics avoid
URL/header credentials. Provider adapters still own their error presentation.

## Validation

`tests/test_progressive_transport.py` runs with standalone Python against a local
threaded HTTP fixture. It covers completion identities, bounded concurrency,
paused backpressure, cancellation, helper exit, retries, gzip modes, truncation
and size limits.

`tests/check_concurrent_tiles.py` runs in background Blender with `--artifacts`
and optionally `--helsinki`. It compares sequential blocking, one/four-worker
blocking and progressive imports, with joining on/off. Nested manifest session
queries, distinct transforms, mixed GLB/B3DM content, out-of-order completion,
partial failures and helper cleanup are checked. It records total wall time,
consumer processing time and peak HTTP concurrency. Delays are controlled local
fixture latency; these measurements are not live-provider performance claims.

`tests/check_progressive_import.py` checks geometry/context/collection/lifecycle
behavior; `tests/check_progressive_ui.py` runs real timers and verifies blocking
and progressive Undo/Redo, cancellation, edit-mode pause and shutdown with the
concurrent downloader. All generated results belong under
`C:/dev/workspace/artifacts/concurrent_downloads`.

### Initial timing (Blender 5.2)

The final local-fixture run used six Helsinki tile instances (the three existing
CC-BY-SA 4.0 fixtures in two transformed branches), with a 350 ms delay for the
first tile and 120 ms for each other tile. These are individual runs, not a
statistical live-network benchmark; helper startup is included.

| Mode, separate objects | Wall time | Peak HTTP requests |
| --- | ---: | ---: |
| Previous sequential blocking path | 1.21 s | 1 |
| Blocking, one helper worker | 1.67 s | 1 |
| Blocking, four helper workers | 0.81 s | 4 |
| Progressive, one helper worker | 1.69 s | 1 |
| Progressive, four helper workers | 0.80 s | 4 |

The progressive runner here was stepped in background Blender for timing; UI
responsiveness and undo were checked separately in real UI processes. Raw timing,
import order, geometry comparison and error results are in
`artifacts/concurrent_downloads/integration_final_5_2/concurrent_checks.json`.
