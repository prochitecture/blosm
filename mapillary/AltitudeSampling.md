# Worldwide Mapillary computed-altitude samples

Sampled on 2026-11-12. All altitude statistics below are in metres and refer to
provider `computed_altitude`, independently of Blender camera Z or our assumed
2 m ground-placement policy. Middle 90% means interpolated 5th–95th percentiles.

| Area | Images | Median | Middle 90% | Full range | Negative values |
|---|---:|---:|---:|---:|---:|
| Karl Marx Allee | 750 | 1.98 | -0.57 to 3.81 | -16.88 to 37.42 | 46 |
| Helsinki | 546 | 1.74 | 0.15 to 3.41 | -1.45 to 6.78 | 22 |
| San Francisco, USA | 500 | 101.20 | 86.87 to 106.77 | 73.28 to 112.59 | 0 |
| La Paz, Bolivia | 500 | 2.18 | -1.10 to 6.22 | -10.14 to 34.97 | 42 |
| Sao Paulo, Brazil | 500 | 2.00 | 0.39 to 3.87 | -14.09 to 8.69 | 21 |
| Cape Town, South Africa | 373 | 1.91 | 0.27 to 3.02 | -0.09 to 3.71 | 3 |
| Kathmandu, Nepal | 318 | 2.00 | 0.85 to 3.45 | 0.68 to 4.59 | 0 |
| Tokyo, Japan | 500 | 1.27 | -5.35 to 3.49 | -21.74 to 19.01 | 114 |
| Sydney, Australia | 500 | 2.00 | -0.30 to 2.13 | -15.79 to 8.65 | 26 |
| Nairobi, Kenya | 500 | 2.00 | 0.14 to 2.70 | -6.19 to 8.68 | 20 |

## Method and scope

The Berlin and Helsinki rows reuse the full metadata sets in the previously saved
750- and 546-camera live test scenes. The eight additional rows use fresh public
coverage tiles at zoom 14 and Graph API requests for `id`, `computed_altitude`,
and original `altitude`. No photographs were downloaded and no cameras created.

Each new area is a roughly 222 m square centred at the coordinates recorded in
`results.json`. Nairobi was expanded to roughly 1.11 km across because the first
222 m square had only two candidates. Coverage image IDs inside each rectangle
were deduplicated, sorted as exact strings and sampled without replacement,
using Python Random seed 20261112. At most 500 IDs per area were requested;
Cape Town and Kathmandu use all available candidates (373 and 318 respectively).
All requested entities were returned with finite numeric computed altitudes.

These are local image samples, not whole-city or worldwide population estimates.
Nearby frames and multiple sequences are not independent observations. Image
sampling does not weight capture sequences equally. The small-area baseline rows
and capped random new samples have different sampling procedures. Coverage-point
filtering uses the quantized vector-tile positions rather than Graph geometry.

The reusable script is `tests/sample_mapillary_altitudes.py`; `--only Nairobi
--half-span .005` reproduces the expanded Nairobi sampling. HTTP uses the shared
transport with four concurrent requests, 32 MiB per coverage resource and 256 KiB
per metadata resource. Metadata batches contain at most 32 IDs. The token is read
from the excluded local development file and never printed or stored in results.

## Interpretation

Most samples cluster near 2 m, but San Francisco has a median of 101.20 m and a
73.28–112.59 m range. These new observations do not support assuming that the
field universally represents actual height above the ground. They also do not,
by themselves, establish its vertical datum or the cause of the differences.
The earlier two-area ground-relative working hypothesis did not generalize.

No runtime camera-placement behavior was changed for this study. Existing
terrain-plus-estimated-height placement remains in the development branch.

## Evidence

Outputs are under `C:/dev/workspace/artifacts/mapillary_import/world_altitudes/`:
`results.json` contains exact bounds, candidate/request counts and individual
sampled IDs with original/computed altitudes. `summary.json` includes the two
baseline rows; `summary.csv` is the exportable table. Initial run progress is in
`../world_altitudes_run.log`, and expanded Nairobi progress in
`nairobi_expanded.log`. Existing baseline evidence is in
`../altitude_statistics/statistics.json`.
