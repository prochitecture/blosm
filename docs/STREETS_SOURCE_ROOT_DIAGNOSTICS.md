# Diagnosing a circular source-root refusal

`corner_source_support.complete_source_roots()` still returns the existing
boolean result. Its optional `on_refusal` callback runs only when the actual
root-coverage predicate rejects a candidate. Successful candidate arithmetic,
source curves, radius choices, source widths and all numerical allowances remain
unchanged.

The callback receives the original Section and directed End, the exposed root
side, the actual transverse parameter interval and its physical length. A full
root audit identifies `side=None`; a facing-sector audit identifies `side=-1`
for the first original LEFT root and `side=1` for the second original RIGHT root.

The source-cycle acceptance predicate uses `raise_source_root_refusal()` to
preserve this reason as a `CornerContactError`. Radius policy still tries the
same bounded STYLE candidates. A failure of root coverage must not be described
as insufficient finite approach length merely because a generic acceptance
callback returned false. A separate original Route-seam refusal is also reported
as such.

Tests are in `tests/test_python_street_corner_source_support.py`. A finite
three-arm case verifies the failed native half/interval/length and unchanged
input geometry. A successful candidate verifies that an installed callback is
not invoked and the original result and geometry remain unchanged. The
diagnostic, compound-contact and ordinary-singleton suites passed 11 tests with
unchanged imported source hashes.

## Actual Bergamo site 1947

The fresh runtime-90 PLAN advanced beyond the earlier source handoff and
dimensional-recipient failures, then stopped at this original site. The acute
890/end0 to 974/end1 circle has finite regular contacts well inside its actual
31.95 m and 28.74 m approaches. The refusal is a 1.482434 m part of the wider
890 LEFT root half, transverse parameters approximately `[0,0.141184]`. Trying
the permitted radius sequence from 2 m to 0.5 m does not remove that exposure.

This is a real physical source-root condition, far above numerical welding
allowances. It is not resolved by the diagnostic. Existing finite cyclic
support searches, their immutable operands, bound reachability and original
recipient-plane staging are separate qualification investigations. No audit
operand, radius, mouth bound or geometry acceptance condition is relaxed by
this reporting change.

Exact read-only evidence is in `cold30_1947_accept_refusals_report.json` under
`C:/dev/workspace/artifacts/streets/beijing_bergamo_fixes_20261114`. It records
every existing permitted radius and the actual failing predicate. Maintained
tests reconstruct their native inputs from source; development pickles and
test snapshots are not recovery dependencies.
