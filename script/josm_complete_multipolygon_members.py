#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File: josm_complete_members.py

import argparse, sys, time
from typing import List, Dict, Tuple
import xml.etree.ElementTree as ET
import requests
from urllib.parse import quote_plus

JOSM = "http://127.0.0.1:8111"
TIMEOUT = 30

def rc_get(path: str):
    r = requests.get(f"{JOSM}{path}", timeout=TIMEOUT)
    r.raise_for_status()
    return r

def josm_ping() -> dict:
    return rc_get("/version").json()

def josm_complete_relation(rel_id: int, relation_members=True, referrers=True):
    """
    Ask JOSM to add the relation (if missing) and download any missing
    members/referrers into the ACTIVE (current) data layer.
    """
    url = f"/load_object?objects={quote_plus('r'+str(rel_id))}"
    if relation_members: url += "&relation_members=true"
    if referrers:        url += "&referrers=true"
    rc_get(url)

def parse_osm_file(path: str) -> Tuple[List[int], Dict[int, dict]]:
    """Return (relation_ids, tags_by_relation_id) from a local .osm/.osm.xml file."""
    try:
        root = ET.parse(path).getroot()
    except ET.ParseError as e:
        raise RuntimeError(f"Invalid OSM XML in {path}: {e}")
    rel_ids, tag_map = [], {}
    for rel in root.findall("relation"):
        rid = int(rel.attrib["id"])
        tags = {t.attrib["k"]: t.attrib["v"] for t in rel.findall("tag")}
        tag_map[rid] = tags
        rel_ids.append(rid)
    return rel_ids, tag_map

def run(tag: str, osm_file: str, wait_after_each: float, dry_run: bool):
    if "=" not in tag:
        raise ValueError("tag must be 'key=value', e.g. natural=wood")
    key, value = tag.split("=", 1)

    print("JOSM RC:", josm_ping())

    rel_ids, tag_map = parse_osm_file(osm_file)
    to_complete = [rid for rid in rel_ids
                   if tag_map[rid].get("type") == "multipolygon" and tag_map[rid].get(key) == value]

    if not to_complete:
        print(f"No multipolygon relations with {tag} found in {osm_file}.")
        return

    print(f"Found {len(to_complete)} multipolygon relations with {tag}.")
    if dry_run:
        print("Dry-run: would request completion for relation IDs:")
        print(", ".join(map(str, to_complete)))
        return

    print("Requesting downloads into the CURRENT JOSM layer …")
    for rid in to_complete:
        josm_complete_relation(rid, relation_members=True, referrers=True)
        if wait_after_each > 0:
            time.sleep(wait_after_each)
    print("Done.")

def main():
    ap = argparse.ArgumentParser(
        description="Download & add incomplete members of tagged multipolygon relations into the CURRENT JOSM layer."
    )
    ap.add_argument("--tag", required=True, help="OSM tag 'key=value', e.g. natural=wood")
    ap.add_argument("--osm-file", required=True, help="Path to local .osm/.osm.xml saved from JOSM")
    ap.add_argument("--wait-after-each", type=float, default=0.0, help="Seconds to wait after each download request")
    ap.add_argument("--dry-run", action="store_true", help="Print relation IDs and exit")
    args = ap.parse_args()
    try:
        run(args.tag, args.osm_file, args.wait_after_each, args.dry_run)
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr); sys.exit(1)

if __name__ == "__main__":
    main()
