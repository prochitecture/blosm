"""Fail-closed acceptance of one completed native street import record."""

GROUND_QUALIFICATION_POLICY = 'SUCCESSFUL_IMPORT_AND_COMPLETE_GROUND_V2'


def ground_qualification_passed(current):
    """READY reports cannot qualify a failed, cancelled or unfinished import."""
    if current.get('error') or current.get('operator') != ['FINISHED']:
        return False
    return complete_ground_result(current.get('import_result'), current.get('import_reports'))


def complete_ground_result(result, reports):
    """Shared acceptance of the final blocking or progressive session result."""
    if not isinstance(result, dict) or result.get('outcome') != 'success' or result.get('errors') != []:
        return False
    if not isinstance(reports, list) or not reports:
        return False
    for report in reports:
        if not isinstance(report, dict) or report.get('ground_coverage') != 'COMPLETE':
            return False
        surfaces = report.get('surfaces')
        if not isinstance(surfaces, list) or not surfaces:
            return False
        if any(not isinstance(surface, dict) or surface.get('status') != 'READY' for surface in surfaces):
            return False
    return True
