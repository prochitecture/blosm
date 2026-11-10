"""Bounded live regions from the Blosm public 3D Tiles dataset list.

URLs/attribution remain in source; downloaded manifests and tiles stay in memory.
These are test regions, not assertions about provider availability or licenses.
"""
DATASETS = {
    'helsinki': dict(url='https://kartta.hel.fi/3d/mesh/Helsinki_2024/tileset.json',
                     bounds=(24.936, 60.168, 24.944, 60.172), geometric_error=30),
    'aarhus': dict(url='https://sdfe-hosting.virtualcitymap.de/Aarhus/datasource-data/Aarhus_Mesh/tileset.json',
                   bounds=(10.200, 56.155, 10.208, 56.159), geometric_error=30),
    'aurillac': dict(url='https://3d.craig.fr/datasets/Aurillac/3dtiles/tileset.json',
                     bounds=(2.437, 44.924, 2.445, 44.929), geometric_error=30),
    'clermont': dict(url='https://3d.craig.fr/datasets/Clermont/3dtiles/tileset.json',
                     bounds=(3.082, 45.775, 3.090, 45.780), geometric_error=60),
    'strasbourg': dict(url='https://s3.eu-west-2.wasabisys.com/ems-sgct-photomaillage/ODACIT/EMS_PM2022/tileset.json',
                       bounds=(7.744, 48.580, 7.752, 48.585), geometric_error=30),
    'agi': dict(url='https://pelican-public.s3.amazonaws.com/3dtiles/agi-hq/tileset.json',
                bounds=(-75.600, 40.036, -75.594, 40.040), geometric_error=3),
    'san_francisco': dict(url='https://demo3d.aerometrex.com/USA/US0005_San_Francisco/San_Francisco_2cm_SL_split/Scene/Cesium_SF_2cm_mSL.json',
                          bounds=(-122.405, 37.784, -122.399, 37.789), geometric_error=2),
    'helsinki_dense': dict(url='https://kartta.hel.fi/3d/mesh/Helsinki_2024/tileset.json',
                           bounds=(24.936, 60.168, 24.944, 60.172), geometric_error=6),
    # The user supplied this real GLB dataset to exercise unsupported extensions.
    # Its transform places it at these coordinates, regardless of its name.
    'svirnas': dict(url='https://raw.githubusercontent.com/WilliamLiu-1997/3D-Tiles-RendererJS-3DGS-Plugin/refs/heads/main/data/SvirnasAlyt-3dtiles/tileset.json',
                    bounds=(142.8335, -38.5828, 142.835, -38.5816), geometric_error=1,
                    expected_failure='KHR_gaussian_splatting'),
}
