#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File: josm_select_degenerate_inners.py

import argparse, sys
from typing import Dict, List, Tuple, Set
import xml.etree.ElementTree as ET
import requests
from urllib.parse import quote_plus
from collections import defaultdict

JOSM = "http://127.0.0.1:8111"
RC_TIMEOUT = 30
SELECT_BATCH = 400

def rc_get(path: str):
    r = requests.get(f"{JOSM}{path}", timeout=RC_TIMEOUT)
    r.raise_for_status()
    return r

def josm_ping() -> dict:
    return rc_get("/version").json()

def josm_select_ways(way_ids: List[int]):
    """Select ways via /zoom?select=... (batched; dummy bbox)."""
    if not way_ids:
        return
    ids = [f"way{w}" for w in way_ids]
    for i in range(0, len(ids), SELECT_BATCH):
        chunk = quote_plus(",".join(ids[i:i+SELECT_BATCH]))
        rc_get(f"/zoom?left=0&right=0&top=0&bottom=0&select={chunk}")

def parse_osm_file(path: str) -> Tuple[Dict[int, List[int]], List[dict]]:
    """Return (ways_map, relations) from a local .osm/.osm.xml file."""
    try:
        root = ET.parse(path).getroot()
    except ET.ParseError as e:
        raise RuntimeError(f"Invalid OSM XML in {path}: {e}")

    ways: Dict[int, List[int]] = {}
    for w in root.findall("way"):
        wid = int(w.attrib["id"])
        nds = [int(nd.attrib["ref"]) for nd in w.findall("nd")]
        ways[wid] = nds

    relations: List[dict] = []
    for rel in root.findall("relation"):
        tags = {t.attrib["k"]: t.attrib["v"] for t in rel.findall("tag")}
        members = []
        for m in rel.findall("member"):
            if m.attrib.get("type") == "way":
                members.append(("way", int(m.attrib["ref"]), m.attrib.get("role","")))
        relations.append({"id": int(rel.attrib["id"]), "tags": tags, "members": members})
    return ways, relations

def touching_inners(rel: dict, ways: Dict[int, List[int]]) -> Set[int]:
    """
    Mark inner ways as 'degenerate' if:
      • they share ≥1 node with another inner way, OR
      • they share ≥1 node with any outer way.
    Uses a node->frequency map to correctly detect sharing with *other* inners.
    """
    inners = [ref for t,ref,role in rel["members"] if t=="way" and role=="inner" and ref in ways]
    outers = [ref for t,ref,role in rel["members"] if t=="way" and role=="outer" and ref in ways]

    # Count nodes across inner ways (dedup nodes within each way)
    node_counts = defaultdict(int)
    inner_nodesets: Dict[int, Set[int]] = {}
    for wid in inners:
        ns = set(ways[wid])
        inner_nodesets[wid] = ns
        for n in ns:
            node_counts[n] += 1

    outer_nodes = set()
    for ow in outers:
        outer_nodes.update(ways[ow])

    bad: Set[int] = set()
    for wid, ns in inner_nodesets.items():
        # shared with another inner (freq>1 means some other inner also has this node)
        if any(node_counts[n] > 1 for n in ns):
            bad.add(wid)
            continue
        # or touching the outer
        if ns & outer_nodes:
            bad.add(wid)
    return bad

def run(tag: str, osm_file: str):
    if "=" not in tag:
        raise ValueError("tag must be 'key=value', e.g. natural=wood")
    key, value = tag.split("=", 1)

    print("JOSM RC:", josm_ping())

    ways, rels = parse_osm_file(osm_file)
    targets = [r for r in rels if r["tags"].get("type") == "multipolygon" and r["tags"].get(key) == value]
    if not targets:
        print(f"No multipolygon relations with {tag} found in {osm_file}.")
        return

    print(f"Found {len(targets)} multipolygon relations with {tag}.")
    offenders: Set[int] = set()
    for r in targets:
        offenders |= touching_inners(r, ways)

    if offenders:
        arr = sorted(offenders)
        print(f"Selecting {len(arr)} inner rings in JOSM (only objects present in the active layer will highlight)…")
        josm_select_ways(arr)
        print("Done.")
    else:
        print("No touching/degenerate inner polygons found in this file.")

def main():
    ap = argparse.ArgumentParser(
        description="Detect & select degenerate inner polygons from a local .osm file (no downloads)."
    )
    ap.add_argument("--tag", required=True, help="OSM tag 'key=value', e.g. natural=wood")
    ap.add_argument("--osm-file", required=True, help="Path to local .osm/.osm.xml saved from JOSM")
    args = ap.parse_args()
    try:
        run(args.tag, args.osm_file)
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr); sys.exit(1)

if __name__ == "__main__":
    main()
