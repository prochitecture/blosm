"""Read actual rounded-cap parity before constructing Bundle entrance strips."""
from dataclasses import dataclass
import numpy as np
from .bundle_axis import AxialChartUnavailable
from .bundles import BundleProvider, Member
from .model import End, Key
from .providers import CarriagewayProvider, RoadSampling
from .routes import RouteTopology


@dataclass
class SourceRimSampling:
    """The carriageway's shared span/sampling policy, without building its mesh."""
    attachments: dict

    def __call__(self, context, section_id, side):
        section = context.network.section_index[section_id]
        path = context[Key('path', section_id)]
        metadata = CarriagewayProvider.span_metadata(context, section, self.attachments[section_id])
        if metadata['end'] <= metadata['start']+1e-7 or metadata['transition_spans']:
            raise AxialChartUnavailable('Angular cap sampling needs an exposed connecting crossroad span')
        station = RoadSampling.from_metadata(metadata).samples(path, context.network.settings)
        xy = path.point(station, section.offset+side*.5*section.width)
        profile = context[Key('profile', section_id)]
        z = profile.datum+context.network.settings.exaggeration*float(profile.value(side))
        return np.column_stack((xy, np.full(len(xy), z)))


def cap_parities(context, bundle):
    members = []
    for member in bundle.members:
        first, last = member.path[0], member.path[-1]
        ends = (End(first.section_id, 0 if first.forward else 1), End(last.section_id, 1 if last.forward else 0))
        members.append(Member((), np.empty((0, 3)), ends))
    sampling = SourceRimSampling(RouteTopology(context.network).attachments(context.network))
    caps = [BundleProvider.cap(context, *members, end, rim_source=sampling) for end in (0, 1)]
    if any(cap is None for cap in caps):
        raise AxialChartUnavailable('Angular sampling requires two actual junction-facing separator caps')
    return tuple((len(cap.rim)-1)%2 for cap in caps)
