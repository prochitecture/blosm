"""Inspect real OSM preparation before the Python mesh backend, without repairs.

Run in an isolated factory Blender using -- --artifacts DIR --osm FILE [...].
The normal parser/generator and adapters run; numerical task declarations are
checked but mesh construction/publication is intentionally not performed.
This is evidence collection, not a claim that observed geometry is qualified.
"""
import argparse
from collections import Counter, defaultdict
from contextlib import nullcontext
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import pickle
import sys
import time
import traceback
from unittest.mock import patch

import bpy

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_street_progressive as support
from street_issue_qualification import GROUND_QUALIFICATION_POLICY, ground_qualification_passed


def identity(item):
    return dict(kind=type(item).__name__, id=getattr(item, 'id', None))


def inspect(manager, renderer, result):
    active = tuple(manager.iterStreets())
    allowed = {id(s) for s in active} | {id(s) for s in manager.physicalStreets.values()}
    findings, chain_items, sections = [], {}, defaultdict(list)

    def finding(code, **details):
        findings.append(dict(code=code, **details))

    for street in active:
        items, seen, item = [], set(), street.head
        while item is not None and id(item) not in seen:
            seen.add(id(item))
            items.append(item)
            if getattr(item, 'street', None) is not street:
                finding('ITEM_OWNER', street=street.id, item=identity(item),
                        owner=identity(getattr(item, 'street', None)))
            if type(item).__name__ == 'Section':
                sections[id(item)].append(street.id)
                for name in ('width', 'offset'):
                    if not isinstance(getattr(item, name, None), (int, float)):
                        finding('DIMENSION_TYPE', street=street.id, item=identity(item), field=name)
            successor = getattr(item, 'succ', None)
            if successor is not None and getattr(successor, 'pred', None) is not item:
                finding('CHAIN_RECIPROCAL', street=street.id, item=identity(item),
                    successor=identity(successor), successor_pred=identity(getattr(successor, 'pred', None)))
            item = successor
        if item is not None:
            finding('CHAIN_CYCLE', street=street.id, at=identity(item))
        if not items or items[-1] is not street.tail:
            finding('CHAIN_TAIL', street=street.id)
        if street.head is not None and getattr(street.head, 'pred', None) is not None:
            finding('CHAIN_HEAD', street=street.id)
        chain_items[street.id] = items
        for end in ('pred', 'succ'):
            connector = getattr(street, end)
            if connector is not None and type(connector).__name__ != 'IntConnector':
                finding('ENDPOINT_TYPE', street=street.id, end=end, value=identity(connector))
    for owners in sections.values():
        if len(set(owners)) > 1:
            finding('DUPLICATE_SECTION', streets=owners)

    rings = 0
    for intersection in manager.majorIntersections.values():
        first, connector, seen = intersection.startConnector, intersection.startConnector, set()
        if first is None:
            continue
        rings += 1
        while connector is not None and id(connector) not in seen:
            seen.add(id(connector))
            successor = getattr(connector, 'succ', None)
            if successor is None or successor.pred is not connector:
                finding('RING_RECIPROCAL', intersection=intersection.id, connector=connector.id)
            if connector.intersection is not intersection:
                finding('RING_OWNER', intersection=intersection.id, connector=connector.id)
            item = connector.item
            if type(item).__name__ == 'Street' and id(item) not in allowed:
                finding('GHOST_CONNECTOR_STREET', intersection=intersection.id, street=item.id)
            connector = successor
        if connector is not first:
            finding('RING_CYCLE', intersection=intersection.id)

    missing_members = []
    for bundle in manager.iterBundles():
        for member in set(bundle.streetsHead + bundle.streetsTail):
            if id(member) not in allowed:
                finding('BUNDLE_MEMBER_UNAVAILABLE', bundle=bundle.id, street=member.id)
            elif member.id not in manager.streets:
                missing_members.append(dict(bundle=bundle.id, street=member.id, retained=True))

    propagation_checks = 0
    for transition in (*manager.transitionSideLanes, *manager.transitionSymLanes):
        item = transition.succ
        if item is None:
            finding('TRANSITION_NO_SUCCESSOR', item=identity(transition))
            continue
        width, offset, seen = item.width, item.offset, set()
        item = item.succ
        while item is not None and id(item) not in seen:
            seen.add(id(item))
            if type(item).__name__ == 'Section':
                nominal = getattr(item, 'sourceWidth', None)
                nominal = item.width if nominal is None else nominal
                reference = transition.succ
                original = getattr(reference, 'sourceWidth', None)
                original = reference.width if original is None else original
                if (abs(nominal-original) > 1e-7
                        or item.totalLanes != reference.totalLanes):
                    break  # A real original dimensional boundary stops correction propagation.
                propagation_checks += 1
                wrong = abs(item.width - width) > 1e-7
                if type(transition).__name__ == 'SideLane':
                    wrong |= abs(item.offset - offset) > 1e-7
                if wrong:
                    finding('TRANSITION_PROPAGATION', transition=identity(transition),
                        item=identity(item), expected=[width, offset], actual=[item.width, item.offset])
            elif not (type(item).__name__ == 'Intersection' and item.isMinor):
                break
            item = item.succ

    # Record source-backed examples, not transient numeric IDs alone.
    for row in findings:
        sid = row.get('street')
        if sid in chain_items:
            row['osm_way_ids'] = sorted({way for item in chain_items[sid]
                                        for way in getattr(item, 'sourceWayIds', ())})
    result.update(dict(streets=len(active), sections=len(sections), bundles=len(manager.bundles),
        major_rings=rings, right_hand_traffic=manager.rightHandTraffic,
        transition_propagation_checks=propagation_checks,
        findings=findings, finding_counts=dict(Counter(row['code'] for row in findings)),
        retained_bundle_members=missing_members,
        xy_audit=getattr(manager, 'streetPreparationAudit', {})))

    snapshot = renderer.adapter.read(manager)
    result['snapshot_issues'] = snapshot.issues
    result['ownership_audit'] = snapshot.ownership_audit
    from blosm.way.street_network.python_backend.adapter import components
    from blosm.way.street_network.physical_sources import PhysicalStreetInventory
    from blosm.way.street_network.python_backend.model import Settings
    from blosm.way.street_network.python_backend.backend import PythonMeshBackend
    from blosm.way.street_network.python_backend.updates import DependencyIndex
    settings = Settings(renderer.design.cell_size, renderer.design.resolution, renderer.design.exaggeration)
    physical = PhysicalStreetInventory.view(snapshot)
    adapted = renderer.converter.read(physical, 'audit', settings,
        rise=renderer.design.profile_rise, interpolation=renderer.design.interpolation)
    if getattr(support.ARGS, 'export_inputs', False):
        # Preserve inputs for small offline replays; these are generated evidence,
        # not repository fixtures and never replace the original OSM source.
        target = support.ARGS.artifacts / Path(renderer.app.osmFilepath).stem
        target.mkdir(exist_ok=True)
        (target / 'snapshot.json').write_text(json.dumps(asdict(snapshot), indent=2), encoding='utf-8')
        (target / 'network.json').write_text(json.dumps(adapted.network.to_dict(), indent=2), encoding='utf-8')
        result['input_capture'] = str(target)
    result['adapter_issues'] = adapted.issues
    result['components'] = []
    supported_inputs = []
    source_support = None
    if getattr(support.ARGS, 'source_support', False):
        from blosm.way.street_network.python_backend.source_support import ClippedApproachSupport
        source_support = ClippedApproachSupport(physical, snapshot.clipped_source_context)
    source_ways = {item['id']: item.get('source_way_ids', ()) for street in snapshot.streets.values()
                   for item in street['items'] if item['kind'] == 'Section'}
    for network in components(adapted.network):
        ids = {s.id for s in network.sections}
        row = dict(sections=len(ids), blocked=bool(ids & adapted.blocked_sections))
        if source_support is not None:
            preparation = source_support.prepare(network)
            row['source_support'] = preparation.report
            row['source_support_error'] = preparation.error
            network = preparation.network
            supported_inputs.append(network.to_dict())
        if not row['blocked']:
            try:
                engine = PythonMeshBackend(network.settings)
                tasks = tuple(task for provider in engine.providers for task in provider.tasks(network))
                row['repeated_dependencies'] = [dict(task=task.key.text(),
                    repeated=[key.text() for key, count in Counter(task.dependencies).items() if count > 1])
                    for task in tasks if len(set(task.dependencies)) != len(task.dependencies)]
                DependencyIndex(tasks)
                row['tasks'] = len(tasks)
            except Exception as exc:
                row['task_error'] = str(exc)
                row['osm_way_ids'] = sorted({way for sid in ids for way in source_ways.get(sid, ())})
        result['components'].append(row)
    if supported_inputs and getattr(support.ARGS, 'export_inputs', False):
        (target / 'supported_networks.json').write_text(json.dumps(supported_inputs, indent=2), encoding='utf-8')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--artifacts', required=True, type=Path)
    parser.add_argument('--osm', required=True, action='append', type=Path)
    parser.add_argument('--railways', action='store_true', help='Include rail/tram ways in shared preparation')
    parser.add_argument('--full-mesh', action='store_true', help='Run unmodified construction instead of preparation-only inspection')
    parser.add_argument('--export-inputs', action='store_true', help='Save prepared inputs for offline diagnostics')
    parser.add_argument('--source-support', action='store_true',
        help='Also run normal original-source boundary support before offline construction diagnostics')
    parser.add_argument('--performance-request', action='store_true',
        help='Enable the normal temporary import performance preference in the isolated process')
    parser.add_argument('--require-ground-coverage', action='store_true',
        help='Fail qualification unless every retained ground component is published')
    parser.add_argument('--save-scene', action='store_true', help='Save generated inspection scenes under artifacts')
    parser.add_argument('--capture-failed-products', action='store_true',
        help='Save actual pending numerical products and digests after construction fails, for diagnostics')
    args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:])
    if (args.require_ground_coverage or args.save_scene) and not args.full_mesh:
        parser.error('--require-ground-coverage and --save-scene require --full-mesh')
    args.artifacts.mkdir(parents=True, exist_ok=True)
    (args.artifacts / 'data').mkdir(exist_ok=True)
    support.ARGS = args
    support.FIXTURE = args.osm[0]
    support.enable()
    from blosm.app import osm_import
    from blosm.progressive import blender, ownership
    from blosm.way.street_network.python_backend import scene, renderer
    from blosm.way.street_network.python_backend import VERSION
    from blosm.way.street_network.python_backend.updates import BuildContext, UpdateCoordinator
    support.runtime, support.driver = osm_import, blender
    support.scene_runtime, support.ownership = scene, ownership
    checks = support.StreetProgressiveChecks()
    results = []
    current = None

    original_prepare_steps = UpdateCoordinator.prepare_steps

    def capture_failed_products(self, *positional, **keywords):
        try:
            return (yield from original_prepare_steps(self, *positional, **keywords))
        except Exception as error:
            # This observes the normal failed evaluation. Never rebuild digests,
            # substitute products, commit the pending update or suppress failure.
            trace = error.__traceback__
            context = None
            while trace is not None:
                candidate = trace.tb_frame.f_locals.get('context')
                if isinstance(candidate, BuildContext):
                    context = candidate
                    break
                trace = trace.tb_next
            if context is not None:
                folder = args.artifacts / Path(current['osm']).stem
                folder.mkdir(exist_ok=True)
                target = folder / ('failed_products_%03d.pkl' %
                    len(current.setdefault('failed_product_captures', [])))
                try:
                    with target.open('wb') as stream:
                        pickle.dump(dict(implementation=VERSION, network=context.network,
                            products=context.products, task_digests=context.task_digests,
                            source_digests=context.source_digests,
                            validated_product_keys=tuple(context.task_digests),error=str(error)), stream)
                    current['failed_product_captures'].append(str(target))
                    print('ISSUE_AUDIT_FAILED_PRODUCTS', str(target), flush=True)
                except Exception as capture_error:
                    current.setdefault('capture_errors', []).append(str(capture_error))
            raise

    def audit_render(self, manager, data, session=None):
        current['audit'] = {}
        inspect(manager, self, current['audit'])
        yield None

    replacement = nullcontext() if args.full_mesh else patch.object(
        renderer.PythonStreetRenderer, '_render_steps', audit_render)
    capture = patch.object(UpdateCoordinator, 'prepare_steps', capture_failed_products) if (
        args.full_mesh and args.capture_failed_products) else nullcontext()
    with support.offline(), replacement, capture:
        for path in args.osm:
            start = time.perf_counter()
            current = dict(osm=str(path.resolve()), sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                           railways=args.railways, full_mesh=args.full_mesh,
                           performance_request=args.performance_request)
            print('ISSUE_AUDIT_START', str(path), flush=True)
            if args.full_mesh:
                # A failure before a new session finishes must not reuse the
                # preceding input's successful process-global result.
                osm_import.last_result = None
            try:
                settings = checks.prepare(False, path)
                settings.railways = args.railways
                bpy.context.preferences.addons['blosm'].preferences.preferPerformanceDuringImport = args.performance_request
                current['operator'] = sorted(bpy.ops.blosm.import_data())
                if not args.full_mesh and 'audit' not in current:
                    raise RuntimeError('Importer did not reach the audit renderer')
            except Exception:
                current['error'] = traceback.format_exc()
            if args.full_mesh:
                current['import_result'] = dict(osm_import.last_result or {})
                names = {c.get('blosm_street_report') for c in bpy.data.collections if c.get('blosm_street_report')}
                current['import_reports'] = [json.loads(bpy.data.texts[name].as_string()) for name in names]
                if args.export_inputs:
                    target = args.artifacts / path.stem
                    target.mkdir(exist_ok=True)
                    (target / 'numerical_inputs.json').write_text(json.dumps([
                        scene.SourceObjects.read(controller).to_dict() for controller in support.controllers()],
                        indent=2), encoding='utf-8')
                if args.require_ground_coverage:
                    current['ground_qualification_policy'] = GROUND_QUALIFICATION_POLICY
                    current['ground_qualification_passed'] = ground_qualification_passed(current)
                if args.save_scene:
                    bpy.context.preferences.filepaths.file_preview_type = 'NONE'
                    current['scene_qualified'] = bool(current.get('ground_qualification_passed'))
                    suffix = '' if current['scene_qualified'] else '_unqualified'
                    target = args.artifacts / (path.stem + suffix + '.blend')
                    bpy.ops.wm.save_as_mainfile(filepath=str(target))
                    current['scene'] = str(target)
            current['seconds'] = time.perf_counter() - start
            results.append(current)
            (args.artifacts / 'preparation_audit.json').write_text(json.dumps(results, indent=2), encoding='utf-8')
            audit = current.get('audit', {})
            print('ISSUE_AUDIT_RESULT', json.dumps(dict(osm=str(path), seconds=current['seconds'],
                error=current.get('error'), findings=audit.get('finding_counts'),
                xy_counts=audit.get('xy_audit', {}).get('counts'),
                components=len(audit.get('components', [])),
                outcome=current.get('import_result', {}).get('outcome'))), flush=True)
    if args.require_ground_coverage and not all(r.get('ground_qualification_passed') for r in results):
        raise AssertionError('Ground street coverage qualification failed; see preparation_audit.json')


if __name__ == '__main__':
    main()
