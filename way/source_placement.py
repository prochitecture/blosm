"""Interpret an OSM way reference without moving its source controls.

Offsets are positive to the right of the prepared Section. Constant one-way
placement is supported; a variable or ambiguous declaration remains explicit
unresolved input rather than silently acquiring a centered carriageway.
"""
import math
import re


def placement_signature(tags):
    return tuple(sorted((key, value) for key, value in tags.items()
                        if key == 'placement' or key.startswith('placement:')))


def constant_placement(tags, width, lanes, *, direction=1):
    signature = placement_signature(tags)
    if not signature:
        return dict(status='NOT_TAGGED', offset=0.)

    def unresolved(reason):
        return dict(status='UNRESOLVED', reason=reason, tags=signature)

    if (type(direction) is not int or direction not in (-1, 1) or type(lanes) is not int or lanes < 1
            or not math.isfinite(width) or width <= 0):
        return unresolved('Placement requires finite dimensions and one original way direction')
    if tags.get('oneway') not in ('yes', 'true', '1', '-1'):
        return unresolved('Constant placement currently requires an explicit one-way road')
    if any(key not in ('placement', 'placement:start', 'placement:end')
           for key, _ in signature):
        return unresolved('Directional placement needs a separate directional lane model')
    base = tags.get('placement')
    values = (tags.get('placement:start', base), tags.get('placement:end', base))
    if values[0] is None or values[0] != values[1] or (base is not None and base != values[0]):
        return unresolved('Variable placement needs an explicit dimensional transition')
    value = values[0]
    if not isinstance(value, str):
        return unresolved('Constant placement must be an OSM tag value')
    match = re.fullmatch(r'(left_of|middle_of|right_of):([1-9][0-9]*)', value)
    if match is None:
        return unresolved('Unsupported constant placement value')
    lane = int(match[2])
    if lane > lanes:
        return unresolved('Placement lane is outside the resolved carriageway')

    widths = [width / lanes] * lanes
    if 'width:lanes' in tags:
        fields = str(tags['width:lanes']).split('|')
        if (len(fields) != lanes or any(re.fullmatch(r'[0-9]+(?:\.[0-9]+)?(?: m)?', v) is None
                                      for v in fields)):
            return unresolved('Placement needs complete metric lane widths')
        widths = [float(v.removesuffix(' m')) for v in fields]
        if (any(v <= 0 for v in widths) or not math.isclose(sum(widths), width, abs_tol=1e-6, rel_tol=1e-8)):
            return unresolved('Placement lane widths disagree with the resolved carriageway width')

    # Lane numbers increase from left to right in the original OSM way frame.
    # First locate the way within that full width, then locate the road center
    # relative to the way. Reverse only when the prepared traversal is reversed.
    fraction = {'left_of': 0., 'middle_of': .5, 'right_of': 1.}[match[1]]
    reference = sum(widths[:lane-1]) + fraction * widths[lane-1]
    return dict(status='RESOLVED', offset=direction * (.5 * width - reference),
                tags=signature, original_direction=direction, lane_widths=tuple(widths),
                reference_from_left=reference, policy='CONSTANT_ONEWAY_PLACEMENT_V1')
