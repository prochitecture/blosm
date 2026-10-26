"""Authored inputs for the complete Python street inspection gallery.

Definitions use the shared model and fixture factories, without Blender scene
state or saved output meshes. The legacy catalogue stays available for existing
regressions; networks() selects the canonical combined-gallery variants.
"""
from dataclasses import replace
from _python_angular_fixtures import angular_fixture_networks, angular_chain_fixture_networks
from _python_street_fixtures import (fixture_networks, ownership_fixture_networks, height_fixture_networks,
                                    transition_fixture_networks, clearance_fixture_networks, side_junction_fixture_networks,
                                    interior_fixture_networks, consumed_cap_fixture_networks, locality_fixture_networks,
                                    terminal_fixture_networks, radius_fixture_networks, source_support_fixture_networks,
                                    attached_overlap_fixture_networks, multisection_attached_fixture_networks,
                                    skew_bundle_fixture_networks, dimensional_bundle_fixture_networks,
                                    one_sided_free_cap_fixture_networks, absorbed_transition_fixture_networks,
                                    nonparallel_terminal_fixture, higher_valence_fixture_networks, absorbed_profile_fixture_networks,
                                    two_merged_terminal_fixture, exact_bent_seam_fixture, sharp_exact_seam_fixture,
                                    compound_seam_fixture_networks, compound_internal_crown_fixture,
                                    compound_route_fixture_networks, sharp_smooth_seam_fixture)
import _python_angular_fixtures as angular
from _python_mixed_fixtures import (mixed_transition_networks,translated_corner_seed_network,
                                   uniform_mixed_networks,uniform_mixed_bundle_network,sharp_uniform_mixed_networks)

def base_networks(model):
    """Original core fixture families, before the angular extensions."""
    return (*fixture_networks(model), *ownership_fixture_networks(model),
                                *height_fixture_networks(model), *transition_fixture_networks(model),
                                *clearance_fixture_networks(model), *side_junction_fixture_networks(model),
                                *interior_fixture_networks(model), *consumed_cap_fixture_networks(model),
                                *locality_fixture_networks(model), *terminal_fixture_networks(model),
                                *radius_fixture_networks(model), *source_support_fixture_networks(model),
                                *attached_overlap_fixture_networks(model), *multisection_attached_fixture_networks(model),
                                *skew_bundle_fixture_networks(model), *dimensional_bundle_fixture_networks(model),
                                *one_sided_free_cap_fixture_networks(model), *absorbed_transition_fixture_networks(model),
                                nonparallel_terminal_fixture(model), *higher_valence_fixture_networks(model),
                                two_merged_terminal_fixture(model), two_merged_terminal_fixture(model,bent_seams=True),
                                exact_bent_seam_fixture(model), sharp_exact_seam_fixture(model,length=30.),
                                sharp_exact_seam_fixture(model,length=60.), *absorbed_profile_fixture_networks(model),
                                *compound_seam_fixture_networks(model), compound_internal_crown_fixture(model),
                                *compound_route_fixture_networks(model), sharp_smooth_seam_fixture(model))

def legacy_networks(model):
    """Preserve the original 86 raw fixture definitions for comparison."""
    return (*base_networks(model),*angular_fixture_networks(model),*angular_chain_fixture_networks(model))

def additional_networks(model):
    """Later island, attached and internal-crossing variants, with anchors."""
    return (
        *(angular.angular_island_chain_network(model,mode) for mode in (0,3)),
        *(angular.attached_angular_network(model,mode) for mode in (0,3)),
        *(angular.attached_angular_state_network(model,mode,state=state)
          for state in ('separated','multiple') for mode in (0,3)),
        *(angular.attached_angular_network(model,mode,radius=20.) for mode in (0,3)),
        *(replace(source,furniture=(*source.furniture,model.Placement(904,301,5.)))
          for state in ('separated','multiple') for mode in (0,3)
          for source in (angular.internal_angular_state_network(model,mode,state=state),)),
        *(replace(source,furniture=(*source.furniture,model.Placement(904,301,2.,offset=-2.5)))
          for mode in (0,3)
          for source in (angular.internal_angular_state_network(model,mode,state='overlap'),)),
    )

def networks(model):
    """Complete gallery, including dimensional seams in the curved chains."""
    # Match the canonical saved gallery: reversed middle Sections are 7 m wide.
    # The raw legacy factories retain their original 6 m widths.
    original=legacy_networks(model)
    prepared=[]
    for network in original:
        if network.uuid.startswith('angular_chain_'):
            network=replace(network,sections=tuple(replace(section,width=7.)
                if section.id in (117,229) else section for section in network.sections))
        prepared.append(network)
    return (*prepared,*additional_networks(model),translated_corner_seed_network(model),
            *mixed_transition_networks(model),*uniform_mixed_networks(model),uniform_mixed_bundle_network(model),
            uniform_mixed_bundle_network(model,overlap=True),*sharp_uniform_mixed_networks(model))
