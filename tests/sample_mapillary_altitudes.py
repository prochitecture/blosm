"""Bounded, metadata-only geographic sampling; no Blender or photo downloads."""
import hashlib
import argparse
import json
import math
from pathlib import Path
import random
import statistics
import sys
import time

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from downloads import Request
from downloads.transport import DownloadTransport
from downloads.diagnostics import redact_text
from progressive.core import Job, Download, Drain
from mapillary.tiles import image_points, tile_coordinates, contains, identifier

OUT = Path('C:/dev/workspace/artifacts/mapillary_import/world_altitudes')
OUT.mkdir(parents=True, exist_ok=True)
token = Path('C:/dev/workspace/blosm/mapillary/mapillary_token.txt').read_text(encoding='utf-8-sig').strip()
assert token and not any(c.isspace() for c in token), 'Invalid development token'
headers = {'Authorization': 'OAuth ' + token}
scope = hashlib.sha256(token.encode()).hexdigest()
CITIES = [('San Francisco, USA', -122.414, 37.793), ('La Paz, Bolivia', -68.134, -16.496),
          ('Sao Paulo, Brazil', -46.655, -23.561), ('Cape Town, South Africa', 18.415, -33.921),
          ('Nairobi, Kenya', 36.821, -1.286), ('Kathmandu, Nepal', 85.307, 27.704),
          ('Tokyo, Japan', 139.701, 35.659), ('Sydney, Australia', 151.209, -33.86)]


def stats(values):
    values = sorted(values)
    if not values: return {'count': 0}
    def q(p):
        x = (len(values)-1)*p
        i = int(x)
        return values[i]+(values[min(i+1,len(values)-1)]-values[i])*(x-i)
    return dict(count=len(values), minimum=values[0], p05=q(.05), median=statistics.median(values),
                p95=q(.95), maximum=values[-1], negative=sum(v<0 for v in values))


def drive(steps, transport):
    job = Job(steps, transport)
    start = time.monotonic()
    while not job.done:
        job.step()
        if time.monotonic()-start > 240:
            job.cancel()
            raise RuntimeError('Sample stage exceeded four-minute deadline')
        if job.waiting: time.sleep(.001)
    if job.error: raise job.error


parser = argparse.ArgumentParser()
parser.add_argument('--only', help='Resample one matching area, preserving the other results')
parser.add_argument('--half-span', type=float, default=.001, help='Latitude half-span in degrees')
args = parser.parse_args()
results = json.loads((OUT/'results.json').read_text())['areas'] if args.only else []
for city, lon, lat in CITIES:
    if args.only and args.only.lower() not in city.lower(): continue
    delta_lon = args.half_span/math.cos(math.radians(lat))
    bounds = (lon-delta_lon, lat-args.half_span, lon+delta_lon, lat+args.half_span)
    row = dict(area=city, bounds=bounds, target_sample=500)
    ids = set()
    def consume_tile(tile, data, error):
        if error: raise error
        for record in image_points(data, tile):
            if record and contains(bounds, *record['coverage_geometry']):
                ids.add(record['id'])
                if len(ids)>100000: raise RuntimeError('Coverage sample ID ceiling exceeded')
            yield None
    def discover():
        for tile in tile_coordinates(bounds):
            z,x,y = tile
            url = f'https://tiles.mapillary.com/maps/vtp/mly1_computed_public/2/{z}/{x}/{y}?access_token={token}'
            yield Download(Request(url, decode_gzip=True, auth_scope=scope),
                           lambda data,error,tile=tile: consume_tile(tile,data,error))
        yield Drain()
    values, originals, samples = [], [], []
    seen = set()
    def consume_metadata(wanted, data, error):
        if error: raise error
        document = json.loads(data)
        if 'error' in document: raise RuntimeError('Graph API returned an error')
        for entry in document.get('data', []):
            identity = identifier(entry.get('id'))
            if identity not in wanted or identity in seen: continue
            seen.add(identity)
            value, original = entry.get('computed_altitude'), entry.get('altitude')
            valid = lambda v: isinstance(v,(int,float)) and not isinstance(v,bool) and math.isfinite(v)
            if valid(value): values.append(value)
            if valid(original): originals.append(original)
            samples.append({'id':identity,'computed_altitude':value,'altitude':original})
            yield None
    try:
        print('Sampling '+city, flush=True)
        drive(discover(), DownloadTransport(max_requests=4,max_bytes=32*1024*1024,max_buffer_bytes=128*1024*1024))
        candidates = sorted(ids)
        selected = random.Random(20261112).sample(candidates,min(500,len(candidates)))
        def metadata():
            for i in range(0,len(selected),32):
                group = selected[i:i+32]
                url = 'https://graph.mapillary.com/images?image_ids='+','.join(group)+'&fields=id,computed_altitude,altitude'
                yield Download(Request(url,headers=headers,auth_scope=scope),
                               lambda data,error,wanted=set(group): consume_metadata(wanted,data,error))
            yield Drain()
        drive(metadata(),DownloadTransport(max_requests=4,max_bytes=256*1024,max_buffer_bytes=1024*1024))
        row.update(status='passed', coverage_candidates=len(ids), requested=len(selected), returned=len(seen),
                   computed_altitude=stats(values),original_altitude=stats(originals),samples=samples)
    except Exception as error:
        row.update(status='failed',error=redact_text(error,headers=headers))
    results = [previous for previous in results if previous['area'] != city] + [row]
    text = redact_text(json.dumps({'date':'2026-11-12','seed':20261112,'areas':results},indent=2),headers=headers)
    assert token not in text, 'Credential found in report'
    (OUT/'results.json').write_text(text,encoding='utf-8')
    print(json.dumps({k:v for k,v in row.items() if k!='samples'}),flush=True)
