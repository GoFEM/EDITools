#!/usr/bin/env python3
"""Download verified release sources; defaults to the two compile-time libraries."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import urllib.error
import urllib.request


def digest(path):
    checksum = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            checksum.update(block)
    return checksum.hexdigest()


def fetch(entry, output):
    target = output / entry['file']
    if target.is_file() and digest(target) == entry['sha256']:
        print(f'Verified {target.name}', flush=True)
        return
    temporary = target.with_suffix(target.suffix + '.part')
    failures = []
    for url in entry['urls']:
        try:
            request = urllib.request.Request(url, headers={'User-Agent': 'EDITools source downloader'})
            with urllib.request.urlopen(request, timeout=60) as source, temporary.open('wb') as destination:
                shutil.copyfileobj(source, destination)
            if digest(temporary) != entry['sha256']:
                raise ValueError('SHA-256 mismatch')
            temporary.replace(target)
            print(f'Downloaded {target.name}', flush=True)
            return
        except (OSError, ValueError, urllib.error.URLError) as error:
            failures.append(f'{url}: {error}')
        finally:
            temporary.unlink(missing_ok=True)
    raise RuntimeError('\n'.join(failures))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('build/sources'))
    parser.add_argument('--manifest', type=Path, default=Path(__file__).resolve().parents[1] / 'packaging/sources-v1.0.json')
    parser.add_argument('--dependencies', action='store_true', help='Also download upstream runtime-library sources (large).')
    args = parser.parse_args()
    catalog = json.loads(args.manifest.read_text())
    entries = catalog['vendor'] + (catalog['dependency_sources'] if args.dependencies else [])
    args.output.mkdir(parents=True, exist_ok=True)
    for entry in entries:
        if Path(entry['file']).name != entry['file']:
            raise ValueError('Source archive name must be a filename')
        fetch(entry, args.output)
