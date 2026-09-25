#!/usr/bin/env python3
"""Pack official Natural Earth 10m SHP ZIPs as Qt-compressed line geometry.

Usage: python3 tools/pack_natural_earth.py /path/to/downloads /path/to/map-data
Downloads retain their official ne_10m_*.zip filenames.
Only geometry is retained; polygons become shoreline rings. No simplification.
Output: little-endian uint32 line count, then per line uint32 point count and
float32 longitude/latitude pairs; prefixed by big-endian uncompressed size and
zlib-compressed (Qt qCompress format). Requires only Python's standard library.
"""
import argparse
import hashlib
import json
from pathlib import Path
import struct
import zipfile
import zlib

LAYERS = {
    'coast': ('physical', 'coastline'),
    'countries': ('cultural', 'admin_0_boundary_lines_land'),
    'states': ('cultural', 'admin_1_states_provinces_lines'),
    'rivers': ('physical', 'rivers_lake_centerlines'),
    'lakes': ('physical', 'lakes'),
}

def pack(source, output):
    with zipfile.ZipFile(source) as archive:
        shp = archive.read(next(n for n in archive.namelist() if n.endswith('.shp')))
        versions = [n for n in archive.namelist() if n.endswith('VERSION.txt')]
        version = archive.read(versions[0]).decode().strip() if versions else 'unspecified'
    lines = []
    offset = 100
    while offset < len(shp):
        size = struct.unpack_from('>I', shp, offset + 4)[0] * 2
        record = shp[offset + 8:offset + 8 + size]
        offset += 8 + size
        kind = struct.unpack_from('<I', record)[0]
        if kind == 0:
            continue
        if kind not in (3, 5):
            raise ValueError(f'Expected polyline or polygon, got {kind}')
        parts, points = struct.unpack_from('<II', record, 36)
        indices = list(struct.unpack_from(f'<{parts}I', record, 44)) + [points]
        start = 44 + parts * 4
        for a, b in zip(indices, indices[1:]):
            coords = struct.unpack_from(f'<{2*(b-a)}d', record, start + a * 16)
            lines.append(struct.pack('<I', b-a) + struct.pack(f'<{len(coords)}f', *coords))
    raw = struct.pack('<I', len(lines)) + b''.join(lines)
    output.write_bytes(struct.pack('>I', len(raw)) + zlib.compress(raw, 9))
    return {'version': version, 'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
            'packed_sha256': hashlib.sha256(output.read_bytes()).hexdigest(), 'lines': len(lines)}

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('downloads', type=Path)
    parser.add_argument('destination', type=Path)
    args = parser.parse_args()
    destination = args.destination
    destination.mkdir(parents=True, exist_ok=True)
    manifest = {}
    for layer, (category, name) in LAYERS.items():
        metadata = pack(args.downloads / f'ne_10m_{name}.zip', destination / f'{layer}.bin')
        metadata['url'] = f'https://naturalearth.s3.amazonaws.com/10m_{category}/ne_10m_{name}.zip'
        manifest[layer] = metadata
    (destination / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
