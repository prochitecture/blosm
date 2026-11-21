# OSM way placement and carriageway offsets

An OSM way does not always represent the middle of a carriageway. The
[`placement` tag](https://wiki.openstreetmap.org/wiki/Proposed_features/placement)
can locate it at a lane center or lane boundary. Blosm retains the original way
controls and uses the Section's signed lateral offset to place the carriageway.
Positive offset means right of the prepared Section direction.

## Constant one-way placement

[`constant_placement`](../way/source_placement.py) interprets `left_of:k`,
`middle_of:k` and `right_of:k` on explicit one-way roads. Lanes are numbered from
left to right in the original OSM way direction. Let the full width be $W$, lane
widths be $w_i$, and the selected lane number be $k$. Set $a$ to zero, one half,
or one for its left edge, middle, or right edge. The way's position measured from
the left carriageway edge is

$$
r = \sum_{i=1}^{k-1} w_i + a w_k.
$$

If $d$ is one for an original traversal and minus one for a reversed traversal,
the carriageway center offset from the unchanged way is

$$
o = d\left(\frac{W}{2} - r\right).
$$

For three 3.5 m lanes, `placement=right_of:1` therefore produces a 1.75 m offset
to the right. It does not move any original source node or change road width.
Equal lane widths follow the current resolved width/lane-count policy. A fully
specified metric `width:lanes` may supply unequal widths, provided their sum
agrees with that full width.

Identical explicit start/end values also describe a constant placement. Variable
placement, `transition`, directional two-way placement, incomplete lane widths,
invalid lane numbers and ambiguous source direction remain explicit unresolved
input. The Python adapter reports them as preparation errors rather than
silently centering the road or treating supported ground input as an intentional
elevated-road exclusion. Automatic inference of variable placement is not yet
implemented.

## Preparation and ownership

Section aggregation retains changes of placement and relevant lane widths,
including original direction reversals. This prevents one merged Section from
losing the reference position of another OSM way. Section source-segment records
determine the prepared/original direction relationship.

Legacy turn-lane offset inference and pixel corrections may not overwrite an
explicit placement declaration. Dimensional propagation stops at a placement or
relevant lane-width boundary. Each Section retains its `sourcePlacement` record
in the prepared snapshot; the numerical backend receives the resulting width
and signed offset through its existing input contract. Changing that offset
changes actual path/profile dependencies. An earlier numerical plan or product
must not be relabeled for the changed source dimensions.

## Evidence and tests

The original Bergamo extract contains four explicitly placed one-way roads;
Beijing contains none. Two affected Bergamo roads previously overlapped by
roughly 0.56 m. Applying their declared placement separates their full analytical
carriageway fibres while keeping original controls and path frames byte-exact.
This is source-preparation evidence, not a complete-city qualification claim.

The focused tests cover lane references, unequal widths, explicit road width,
reversed traversal, unresolved input, aggregation boundaries and protection from
lane-change inference/propagation. Native preparation tests must run in the
chat's isolated factory Blender profile and load this checkout explicitly:

- [Pure placement policy tests](../tests/test_street_source_placement.py).
- [Native placement preparation tests](../tests/test_street_placement_preparation.py).
- [Dimensional propagation tests](../tests/test_street_dimension_propagation.py).
- [Numerical adapter tests](../tests/test_python_street_adapter.py).

Generated diagnostics and current whole-city acceptance status are recorded in
the [Beijing/Bergamo qualification report](STREETS_BEIJING_BERGAMO_FIXES.md).
