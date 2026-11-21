# Binding Bundle spans to published junction mouths

A visible approach index belongs to the layout that produced it. Source-cycle
reconciliation can change the completed physical approach inventory while
preserving every actual street identity. An earlier positional index must not
select an arm in that different inventory.

The real Bergamo failure demonstrates this distinction: the visible registry
identifies original `End(497, 0)` using projected index 2. In the accepted
Product, physical index 1 carries that source mouth; index 2 belongs to source
588. Treating index 2 as persistent identity incorrectly reports a missing
source alias. The correction applies to all source identities, without a city,
Bundle number, coordinate or street-ID special case.

## Identity and binding

| Record | Meaning |
| --- | --- |
| `End(section_id, end)` | Original directed Section endpoint; `end` is 0 or 1 |
| `VisiblePort.approach_index` | Index in the earlier projected source layout |
| `SourceMouth.endpoint` and `original_incidence` | Explicit original incidence or separately identified interior physical cut |
| `(owner, section_id, end)` | Owner-scoped identity of a derived physical cut; it does not create an original End |
| `PhysicalPortBinding.approach_index` | Index resolved within the accepted physical Product's current layout |
| `Product.peaks` / `derived_peaks` | Published source profile stations, in the corresponding identity namespace |

`separator_ports.bind_physical_port()` resolves the actual accepted owner and
selects its physical layout. Ordinary Products without a separate physical
layout use their existing projected layout. It enumerates that small local
approach inventory and requires exactly one explicit `SourceMouth` with the
requested endpoint and the correct original/derived incidence kind.

For a compound Product, the source identity must also have a published profile,
and its native station must match the published station within the existing
`1e-8 m` comparison. Derived bindings additionally retain the exact owner-scoped
alias declaration. Missing profiles, duplicate aliases, nonfinite or changed
stations and a missing derived physical layout are errors. There is no nearest
point search, source relabeling, fallback to an unrelated arm or spatially
inferred connection.

The result carries the accepted Product, owner, current layout, current index,
the matched alias and the complete physical-mouth alias tuple. It is a local
binding, not a new persistent source ID or a replacement for the published
profile samples.

## Consumers and updates

`merged_run_inventory` uses this binding when hydrating original and derived
span ends, resolving a moved owned terminal, and creating the visible physical
mouth. `bundle_terminal_unions` uses the same binding before comparing two
accepted paired-terminal approach indices. Comparing indices is meaningful
only after both bindings identify the same accepted physical owner/layout.

The visible registry retains its existing projected-layout semantics for clients
that read that same projected layout. The correction does not globally reinterpret
an old index, reorder source records or change any generated surface.

`PHYSICAL_PORT_BINDING_POLICY` versions the separator-port and Bundle-plan Tasks.
The latter also projects broad registry dependencies, so its explicit version
must change even when source geometry and local surface stamps stay equal.
Runtime version
`0.1.0.dev92` invalidates earlier implementation fingerprints. Existing local
Bundle signatures retain their published surface digests, so an accepted source
mouth changes through the normal declared dependency graph. No handler, timer,
extra scene object or persistent alias cache is installed.

## Verification and resumption

`tests/test_python_street_physical_port_bindings.py` uses the portable authored
hydrated-run fixture. It covers reordered accepted inventories, an unrelated
old positional index, missing/duplicate aliases, unpublished profiles, changed
native stations, separate original/derived identity spaces and the Task policy
version. Existing separator-port, external-run, assembly and Bundle-plan tests
exercise the complete consumers. Tests do not require a saved development
scene or captured artifact pickle.

The real captured Bundle plan passes an identity-based algorithm diagnostic
without changing its sources or geometry. The maintained candidate-32 normal
Task builder subsequently passes the same local Bundle plan in 73.32 s, with no
replacement builder or injected group, and all source/input hashes unchanged.
Its genuine upstream geometry was restored from candidate 31's validated failure
prefix; this is explicitly a local consumer check, not a fresh coordinator run.

All eight new binding tests and fifteen related separator/external-run/assembly
and mouth-representation methods pass. The expanded 24-method run also includes
one older complete merged-Bundle test whose two subcases fail at boundary-jet
compilation and fixed-height feasibility. Both identical failures reproduce
independently on immutable candidate 31, before the binding change. Subsequent
[paired chord-jet domains](STREETS_SOURCE_OFFSET_RULES.md) and the existing
[direct-quad grass recovery](STREETS_FIXED_GRASS_QUAD_EARS.md) correct those
causes. A scoped numerical candidate then passes all 43 selected binding,
field, complete merged-Bundle, update/restore, query and height tests. Combined
current-source and native qualification remains separate.

A fresh complete coordinator run remains required. See
[the current city qualification report](STREETS_BEIJING_BERGAMO_FIXES.md) for
their actual status; a replay of upstream completed products is not a fresh
complete import.
