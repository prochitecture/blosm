"""Fresh genuine target-Surface dependency closure from native numerical inputs.

The maintained runtime Task must enable its required source-height policy.
This driver never supplies a group, mesh, cut family,
replacement builder or cached fingerprint. All selected Task objects come from
the actual complete runtime registry and are evaluated by UpdateCoordinator.

Launch with this chat's isolated resource-guard host-Python runner. Arguments:
--source IMMUTABLE_ADDON_SNAPSHOT --native-json NATIVE_TRANSPORT_NETWORKS_JSON
--owner 122 --artifacts NEW_ARTIFACT_DIRECTORY [--section-count COUNT].

This is a target dependency-closure witness, not a complete native import,
progressive UI/cancellation test or authorization to publish a partial update.
"""
import argparse
import faulthandler
import hashlib
import json
from pathlib import Path
import pickle
import sys
import time
import traceback


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',required=True,type=Path)
    parser.add_argument('--native-json',required=True,type=Path)
    parser.add_argument('--owner',type=int,default=122)
    parser.add_argument('--section-count',type=int)
    parser.add_argument('--artifacts',required=True,type=Path)
    args=parser.parse_args(sys.argv[sys.argv.index('--')+1:] if '--' in sys.argv else None)
    source=args.source.resolve();out=args.artifacts.resolve()
    artifact_root=Path('C:/dev/workspace/artifacts/streets').resolve()
    if not out.is_relative_to(artifact_root) or out==artifact_root:
        raise ValueError('Target closure reports need a dedicated street artifact directory')
    manifest_path=source.parent/'source_manifest.json'
    manifest=json.loads(manifest_path.read_text())
    if not manifest.get('qualified') or Path(manifest['snapshot']).resolve()!=source:
        raise ValueError('Target closure requires its coherent immutable runtime snapshot')
    def hashes():return {name:hashlib.sha256((source/name).read_bytes()).hexdigest() for name in manifest['files']}
    before=hashes()
    if before!=manifest['files']:raise ValueError('Immutable target runtime differs from its manifest')
    if any(name=='tested_python_streets' or name.startswith('tested_python_streets.') for name in sys.modules):
        raise ValueError('Target qualification must start in a fresh numerical process')
    sys.path.insert(0,str(source/'tests'))
    import _python_street_package
    if not Path(sys.modules['tested_python_streets'].__file__).resolve().is_relative_to(source):
        raise ValueError('Target qualification loaded another numerical implementation')
    import numpy as np
    from tested_python_streets.model import Network,Key,fingerprint
    from tested_python_streets.backend import PythonMeshBackend
    from tested_python_streets.updates import DependencyIndex,UpdateCoordinator,BuildContext,VERSION
    from tested_python_streets.ownership import PLAN,surface_key,CompoundProduct
    from tested_python_streets.compounds import CompoundGroup
    from tested_python_streets.source_body_heights import prepare as body_guides
    from tested_python_streets.publication_slopes import PublicationSlopeGuard
    from tested_python_streets.four_side_charts import quality
    from tested_python_streets.validation import validate_overlaps
    from tested_python_streets.uv import PlanarMapping,apply_mapping,publication_uv_coordinates
    out.mkdir(parents=True,exist_ok=False)
    report=dict(qualified=False,scope='fresh actual target Surface Task and its declared dependency closure',
        whole_native_import_qualified=False,progressive_ui_qualified=False,source=str(source),
        source_manifest_sha256=hashlib.sha256(manifest_path.read_bytes()).hexdigest(),implementation=VERSION,
        source_heights_required=True,restored_products=0,fabricated_fingerprints=False,
        injected_groups=False,injected_meshes=False,replaced_builders=False)
    clock=time.perf_counter();stack=(out/'stacks.log').open('w')
    faulthandler.dump_traceback_later(60.,repeat=True,file=stack)
    def serializable(value):
        if isinstance(value,np.generic):return value.item()
        if isinstance(value,np.ndarray):return value.tolist()
        if isinstance(value,Key):return value.text()
        if isinstance(value,Path):return str(value)
        raise TypeError('Unsupported target report value: '+type(value).__name__)
    def save():
        report['seconds']=time.perf_counter()-clock
        (out/'report.json').write_text(json.dumps(report,indent=2,default=serializable),encoding='utf8')
    def capture(context, *, target_completed, failure_groups=()):
        # Preserve actual runtime evidence before releasing numerical buffers.
        # A failure capture is diagnostic only: it does not qualify the target
        # or authorize injecting unfinished products into a later update.
        completed=tuple(context.task_digests or ())
        state=dict(implementation=VERSION,network=context.network,
            products=dict(context.products),task_digests=dict(context.task_digests or {}),
            source_digests=dict(context.source_digests or {}),
            completed_task_keys=completed,target_key=target,
            target_completed=target_completed,uncommitted=True,
            reusable_as_validated_update=False,
            failure_groups=tuple(failure_groups),
            source_manifest_sha256=report['source_manifest_sha256'])
        path=out/('actual_completed_target.pkl' if target_completed else 'actual_failure_context.pkl')
        temporary=path.with_suffix('.tmp')
        with temporary.open('wb') as stream:pickle.dump(state,stream,pickle.HIGHEST_PROTOCOL)
        temporary.replace(path)
        report['actual_context_capture']=str(path)
    try:
        raw=args.native_json.read_bytes();values=json.loads(raw)
        if isinstance(values,dict):values=values.get('networks',[values])
        candidates=[value for value in values
            if (args.section_count is None or len(value['sections'])==args.section_count)
            and any(site['id']==args.owner for site in value.get('intersections',()))]
        if len(candidates)!=1:raise ValueError('Native inputs do not identify one requested original network')
        network=Network.from_dict(candidates[0]);original=network.to_dict()
        report.update(native_json=str(args.native_json.resolve()),native_json_sha256=hashlib.sha256(raw).hexdigest(),
            canonical_network_fingerprint=fingerprint(network),network_uuid=network.uuid,sections=len(network.sections))
        runtime=PythonMeshBackend(network.settings)
        tasks=tuple(task for provider in runtime.providers for task in provider.tasks(network))
        graph=DependencyIndex(tasks);target=surface_key(Key('intersection',args.owner))
        if target not in graph.tasks:raise ValueError('Current runtime has no requested target Surface Task')
        selected=set();pending=[target]
        while pending:
            key=pending.pop()
            if key in selected or key not in graph.tasks:continue
            selected.add(key);pending.extend(graph.tasks[key].dependencies)
        if PLAN not in selected:raise ValueError('Target closure did not include fresh physical ownership planning')
        chosen=tuple(graph.tasks[key] for key in graph.order if key in selected)
        source_keys=set(network.records())|{Key('settings',0)}
        if any(set(task.dependencies)-selected-source_keys for task in chosen):
            raise ValueError('Selected target closure omits a current declared dependency')
        class ClosureProvider:
            source_kinds=frozenset().union(*(provider.source_kinds for provider in runtime.providers))
            def tasks(self,current):
                if current is not network:raise ValueError('Target closure changed its original input identity')
                return iter(chosen)
        report['task_definitions']=[dict(key=task.key.text(),version=task.version,
            dependencies=[key.text() for key in task.dependencies],
            builder_module=getattr(task.build,'__module__',type(task.build).__module__),
            builder_name=getattr(task.build,'__qualname__',type(task.build).__qualname__),
            signature_module=(getattr(task.signature,'__module__',type(task.signature).__module__)
                if task.signature is not None else None),
            signature_name=(getattr(task.signature,'__qualname__',type(task.signature).__qualname__)
                if task.signature is not None else None)) for task in chosen]
        save();print('FRESH_TARGET_CLOSURE_BEGIN',target.text(),len(chosen),flush=True)
        coordinator=UpdateCoordinator((ClosureProvider(),),runtime.patterns)
        prepared=coordinator.prepare(network,force=True)
        if not {PLAN,target}<=set(prepared.rebuilt):raise ValueError('Target or PLAN was not genuinely built')
        if any(graph.tasks[task.key] is not task for task in chosen):
            raise ValueError('Target closure replaced an actual runtime Task definition')
        resolved=prepared.products[PLAN].owner(Key('intersection',args.owner))
        report['resolved_owner']=resolved.text() if resolved is not None else None
        product=prepared.products[target]
        if not isinstance(product,CompoundProduct):raise ValueError('Requested Task is not the actual compound owner')
        if product.pavement.metadata.get('source_height_conditions') is not True:
            raise ValueError('Actual maintained target Task has not enabled the required source-height policy')
        if network.to_dict()!=original:raise ValueError('Target preparation changed original source inputs')
        context=BuildContext(network,prepared.products,runtime.patterns,prepared.task_fingerprints,prepared.source_fingerprints)
        capture(context,target_completed=True)
        mesh=product.physical_pavement
        edges=np.sort(np.stack((mesh.quads,np.roll(mesh.quads,-1,axis=1)),axis=-1).reshape(-1,2),axis=1)
        edges,counts=np.unique(edges,axis=0,return_counts=True);boundary=np.unique(edges[counts==1])
        guide=body_guides(context,product.group,mesh,boundary)
        if not len(guide.core_nodes):raise ValueError('Target fixture no longer provides a genuine source-body core')
        body_error=float(abs(mesh.xyz[guide.core_nodes,2]-guide.core_targets).max())
        query_error=float(abs(product.field(mesh.xyz[guide.core_nodes,:2])-guide.core_targets).max())
        if max(body_error,query_error)>1e-8:raise ValueError('Target Task lost an original physical body constraint')
        physical=quality(mesh)
        stored=quality(PublicationSlopeGuard.physical_buffer(mesh,product.field.datum,product.exaggeration or 1.))
        if (max(physical['maximum_edge_ratio'],stored['maximum_edge_ratio'])>8.
                or min(physical['minimum_corner_sine'],stored['minimum_corner_sine'])<.03):
            raise ValueError('Target Task fails unchanged physical/stored pavement shape')
        payload={Key('target_pavement',args.owner):mesh};uv=[]
        if product.physical_grass is not None:payload[Key('target_grass',args.owner)]=product.physical_grass
        overlap=validate_overlaps(payload,cell_size=2.)
        for part in (product.pavement,product.grass):
            if part is None:continue
            mapped=apply_mapping(part.mesh,PlanarMapping());values=publication_uv_coordinates(mapped.uv_layers['UVMap'])
            if not np.isfinite(values).all():raise ValueError('Target UV publication contains nonfinite values')
            uv.append(dict(material=part.mesh.material,shape=list(values.shape)))
        points=mesh.xyz[mesh.quads[:min(100,len(mesh.quads))],:2].mean(1)
        expected=product.height(points);product.release_geometry()
        np.testing.assert_array_equal(product.height(points),expected)
        report.update(qualified=True,genuine_target_task=True,fresh_plan=True,task_count=len(chosen),
            rebuilt=[key.text() for key in prepared.rebuilt],source_fingerprints={key.text():value for key,value in prepared.source_fingerprints.items()},
            task_fingerprints={key.text():value for key,value in prepared.task_fingerprints.items()},timings=prepared.timings,
            pavement_quads=len(mesh.quads),grass_quads=sum(len(value.quads) for key,value in payload.items() if key.kind=='target_grass'),
            body_nodes=len(guide.core_nodes),body_error=body_error,physical_query_error=query_error,
            physical_quality=physical,stored_quality=stored,overlap=overlap,uv=uv,
            release_preserves_queries=True,metadata=product.pavement.metadata)
    except Exception as error:
        report.update(qualified=False,error=str(error),traceback=traceback.format_exc());traceback.print_exc()
        seen_errors=set();seen_groups=set();failure_groups=[];failed_context=None;failure=error
        while failure is not None and id(failure) not in seen_errors:
            seen_errors.add(id(failure));trace=failure.__traceback__
            while trace is not None:
                frame=trace.tb_frame
                candidate_context=frame.f_locals.get('context')
                if failed_context is None and isinstance(candidate_context,BuildContext):
                    failed_context=candidate_context
                for name in ('group','original','current','candidate','result'):
                    group=frame.f_locals.get(name)
                    if isinstance(group,CompoundGroup) and id(group) not in seen_groups and len(failure_groups)<16:
                        seen_groups.add(id(group));failure_groups.append(dict(
                            module=frame.f_globals.get('__name__'),function=frame.f_code.co_name,
                            line=trace.tb_lineno,local_name=name,group=group))
                trace=trace.tb_next
            failure=failure.__cause__ or failure.__context__
        if failed_context is not None:
            try:capture(failed_context,target_completed=False,failure_groups=failure_groups)
            except Exception as capture_error:
                report['failure_group_capture_error']=str(capture_error)
                try:capture(failed_context,target_completed=False)
                except Exception as context_error:report['context_capture_error']=str(context_error)
    finally:
        faulthandler.cancel_dump_traceback_later();stack.close()
        after=hashes();report['source_unchanged']=before==after
        report['source_changed']=[name for name in before if before[name]!=after[name]]
        if not report['source_unchanged']:report['qualified']=False
        save();print('TARGET_TASK_RESULT',report['qualified'],report.get('error'),report['seconds'],flush=True)
    return 0 if report['qualified'] else 1


if __name__=='__main__':raise SystemExit(main())
