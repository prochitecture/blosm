# Shared concurrent downloads

`downloads` contains no Blender imports. `DownloadTransport` runs one ordinary
Python helper with a bounded `ThreadPoolExecutor`; all Blender work stays on the
main thread. The helper is reused throughout an import and accepts requests as
new resources are discovered. Neither requests nor results create dataset files.

Inside Blender, transport creation, submission and polling respect
`bpy.app.online_access`. Disabling online access closes the transport at its next
submission or poll, including the isolated worker and pending requests. Ordinary
Python use remains independent of Blender.

## Interface

- `Request(url, headers={}, decode_gzip=False, verify_tls=True, timeout=90.,
  auth_scope='', auth_error_reasons=())` describes an HTTP
  request. Headers and options are copied into the wire message on submission.
  Raw responses are the default so an intrinsically compressed resource, such as
  an HGT gzip file, can be passed unchanged to its consumer.
  `timeout` is the helper-side resource deadline in seconds (up to 3600).
  The optional authentication policy is explained below.
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

The [Mapillary camera adapter](../mapillary/README.md) also uses this transport
in both progressive and blocking modes, with smaller resource limits and bounded
metadata batches. Its decoding and Blender camera commits remain main-thread
cooperative work; they are not performed by the HTTP worker threads.

Concurrent downloads are enabled by default, independently of progressive import.
The toggle and simultaneous-download limit (default four) are internal scene
parameters, hidden from the UI. They remain accessible in Blender's Python Console:

```python
bpy.context.scene.blosm.concurrentDownloads = True
bpy.context.scene.blosm.downloadWorkers = 4  # 1 through 8
```

Disabling `concurrentDownloads` retains the previous sequential path. Setting
`downloadWorkers` to one exercises the concurrent implementation with one
outstanding request. Changes apply to the next import.

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
to twice, with exponential backoff/jitter and numeric/date Retry-After support.
Incomplete fixed-length and chunked response bodies use this same retry budget.
Each retry starts a complete GET with an empty buffer; partial bytes are discarded.
The original resource deadline still applies across attempts and backoff. Invalid
gzip, size violations and glTF parsing errors are not retried as network failures.
HTTP 429 and 503 responses set a cooldown shared by all workers for that origin
(scheme, hostname and port). Starting with the second 429, each further 429 halves
that origin's concurrency limit, down to one, for the rest of this import. Other
origins have independent limits. Requests already in flight finish normally.
A cooldown longer than a resource's remaining deadline returns an HTTP error
asking the user to retry later; it never causes an early retry.

Each resource has a helper-side deadline, normally 90 seconds, covering all
attempts and cooldown waits. A timer shuts down its network sockets on expiry,
including stalled headers and continuously trickling bodies. Socket operations
also have a timeout of at most 30 seconds. OS DNS resolution cannot be forcibly
interrupted by a Python thread. An independent supervisor tracks unfinished
futures and terminates the entire helper if one remains blocked five seconds
beyond its deadline. Its exit code produces a specific stalled-worker error;
the import stops and retains committed geometry. The supervisor runs separately
from both network workers and IPC, so an unresponsive resolver or blocked socket
write to Blender cannot prevent supervision. Finished futures are excluded.
The deadline ends when the resource finishes in the helper. Waiting for Blender
to receive/import completed bytes does not count against it.

A normal deadline, HTTP error, invalid gzip, exhausted incomplete-response retries
or size violation fails that resource; other tiles continue. Failure reading an HTTP error body preserves
the original HTTP status. The helper sends periodic heartbeats; after 15 seconds
of responsive polling with no incoming data, the caller reports an unresponsive
helper. Long Blender import/tool pauses are excluded. Helper exit, malformed IPC
or a failed helper responsiveness check stops the whole job and closes it. There
is no automatic helper restart. On parent disconnect the helper exits immediately
without waiting for network threads, including threads stuck in OS calls.
On Windows, IPC can report a connection reset before the helper's exit status is
visible. The caller allows up to one second for that status, using nonblocking
polls, so a stalled worker is not misreported as a generic socket failure.

Authentication suppression is opt-in. An adapter supplies an opaque `auth_scope`
identifying its credential/session pair and explicit `auth_error_reasons` from
its provider's JSON error responses. HTTP 401, or HTTP 400/403 with a matching
`reason` or `status`, prevents new attempts for that scope at the same origin.
Unrecognized 400/403 responses remain individual resource errors. Other scopes
remain usable; the transport does not refresh credentials or sessions. The Google
3D Tiles adapter supplies this policy with a hash of its credential/session query.
[Google's error documentation](https://developers.google.com/maps/documentation/tile/error_handling)
lists both credential/session problems and unrelated errors under 400/403, so the
status alone is insufficient for classifying those responses.

HTTP failures retain `HTTPError` compatibility, status, a bounded sanitized body,
`attempts`, `auth_failed` and `download_code`. Other resource failures use
`DownloadError` with `code`, `attempts` and `auth_failed`. Human-readable messages
include the attempt count; zero means the request was rejected before network I/O.
URLs in diagnostics retain the host/path but remove user-info, query values and
fragments. Known credential values and credential fields in provider bodies are
redacted. The existing 3D Tiles console/root error paths use the same redaction,
including the sequential blocking path. Per-tile errors still use console details
and one aggregate import result rather than individual popup reports.

The shared transport verifies TLS by default. The 3D Tiles adapter explicitly
retains its pre-existing certificate-verification bypass for compatibility;
changing that existing provider policy is separate work. Provider adapters still
own their error presentation.

## Validation

`tests/test_progressive_transport.py` runs with standalone Python against a local
threaded HTTP fixture. It covers completion identities, bounded concurrency,
paused backpressure, cancellation, helper exit, retries, gzip modes, truncation
and size limits.

`tests/test_download_resilience.py` adds local fault injection for stalled and
trickling headers/bodies, delayed caller polling, shared cooldowns, repeated 429s,
long Retry-After, credential scope suppression, broken/unusual error bodies,
credential redaction and parent disconnect during active downloads.
It also covers unrecoverable DNS stalls, completed responses waiting in IPC,
Windows reset/exit ordering and bounded retries of incomplete transfers.

`tests/check_concurrent_tiles.py` runs in background Blender with `--artifacts`
and optionally `--helsinki`. It compares sequential blocking, one/four-worker
blocking and progressive imports, with joining on/off. Nested manifest session
queries, distinct transforms, mixed GLB/B3DM content, out-of-order completion,
partial failures and helper cleanup are checked. It records total wall time,
consumer processing time and peak HTTP concurrency. Delays are controlled local
fixture latency; these measurements are not live-provider performance claims.
Additional cases deliberately truncate a manifest and a tile on their first
attempt, including real Helsinki payloads, and compare the recovered geometry.

`tests/check_progressive_import.py` checks geometry/context/collection/lifecycle
behavior; `tests/check_progressive_ui.py` runs real timers and verifies blocking
and progressive Undo/Redo, cancellation, edit-mode pause and shutdown with the
concurrent downloader. All generated results belong under
`C:/dev/workspace/artifacts/concurrent_downloads` or, for the resilience checks,
`C:/dev/workspace/artifacts/download_resilience`.

`tests/check_live_tiles.py` compares real provider imports using one/four workers,
blocking/progressive execution and joined/separate objects. It checks response
hashes, geometry, topology counts, packed textures, tile/error counts and cleanup.
`tests/check_live_tiles_ui.py` checks actual UI timers, incremental tiles, grouped
Undo/Redo, cancellation and aggregate glTF error reporting with eight workers.
Their regions and public URLs are in `tests/live_tile_datasets.py`. Reports belong
under `C:/dev/workspace/artifacts/download_live_validation`; downloaded manifests
and tiles are never written. See [LiveValidation.md](LiveValidation.md) for the
observed results and provider limitations.

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
