#!/usr/bin/env python3
"""Assemble a relocatable Qt/GDAL runtime from a release build and Conda prefix.

Build tools only: Python, patchelf (Linux), pefile (Windows dependency inspection).
The prefix must contain the exact libraries used by the compiler. The matching
repository revision must provide build instructions and dependency source locations.
"""
import argparse
from collections import deque
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tarfile
import zipfile


def run(*args):
    return subprocess.check_output([str(a) for a in args], text=True).strip()


def relocate_prefix(path, prefix):
    """Remove the build prefix while preserving offsets of binary C strings."""
    data = path.read_bytes()
    old = str(prefix).encode()
    if old not in data:
        return
    replacement = b'/opt/editools'
    if b'\0' in data:
        parts = data.split(b'\0')
        data = b'\0'.join(part.replace(old, replacement).ljust(len(part), b'\0') for part in parts)
        path.write_bytes(data)
    else:
        path.write_bytes(data.replace(old, replacement))


def archive_licenses(output):
    """Keep notices and dependency provenance in one archive, preserving existing entries."""
    notices = [output / 'LICENSE'] + sorted(p for p in (output / 'licenses').rglob('*') if p.is_file())
    notices += sorted(p for p in (output / 'share').rglob('*') if p.is_file() and
                      p.name.upper().startswith(('LICENSE', 'LICENCE', 'COPYING', 'COPYRIGHT', 'NOTICE')))
    notices += [output / 'DEPENDENCIES.json', output / 'SYSTEM-DEPENDENCIES.txt']
    notices = [p for p in notices if p.is_file()]
    path = output / 'licenses.zip'
    entries = {}
    if path.exists():
        with zipfile.ZipFile(path) as archive:
            entries = {name: archive.read(name) for name in archive.namelist()}
    temporary = path.with_suffix('.zip.tmp')
    with zipfile.ZipFile(temporary, 'w', zipfile.ZIP_DEFLATED) as archive:
        for notice in notices:
            entries[notice.relative_to(output).as_posix()] = notice.read_bytes()
        for name, data in sorted(entries.items()):
            archive.writestr(name, data)
    with zipfile.ZipFile(temporary) as archive:
        if archive.testzip() is not None:
            raise RuntimeError('License archive integrity check failed.')
    temporary.replace(path)
    for notice in notices:
        notice.unlink()
    if (output / 'licenses').exists():
        shutil.rmtree(output / 'licenses')


def write_launch_note(output, platform, source_ref):
    launch = ('Start EDITools.exe (Windows x64).' if platform == 'windows' else
              'Run ./EDITools (Linux x86-64, glibc >= 2.34, X11/XWayland).')
    (output / 'README.txt').write_text(
        f'EDITools\n\nExtract the entire archive. {launch}\n'
        'Keep the bundled files together. No separate installation is needed.\n'
        'Offline geographic map layers are included.\n\n'
        'The main TF view shows XY/YX and XX/YY columns, induction arrows,\n'
        'and phase-tensor ellipses. XY edits also affect XX; YX edits also affect YY.\n'
        'Undo/Redo restores paired and linked masks.\n\n'
        f'Help: https://github.com/GoFEM/EDITools/blob/{source_ref}/docs/README.md\n'
        'Licenses and dependency details: licenses.zip (no extraction needed).\n'
        f'Source and build instructions: https://github.com/GoFEM/EDITools/tree/{source_ref}#building\n', encoding='utf-8')


def copy_map_data(source, output):
    """Include only the verified geometry and provenance from the map packer."""
    manifest = json.loads((source / 'manifest.json').read_text())
    target = output / 'share/natural-earth'
    target.mkdir(parents=True, exist_ok=True)
    for layer in ('coast', 'countries', 'states', 'rivers', 'lakes'):
        path = source / (layer + '.bin')
        if hashlib.sha256(path.read_bytes()).hexdigest() != manifest[layer]['packed_sha256']:
            raise RuntimeError(f'Natural Earth checksum mismatch: {layer}')
        shutil.copy2(path, target / path.name)
    shutil.copy2(source / 'manifest.json', target / 'manifest.json')
    shutil.copy2(source / 'LICENSE.txt', target / 'LICENSE.txt')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--platform', choices=('linux', 'windows'), required=True)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--prefix', type=Path, required=True)
    parser.add_argument('--cache', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--repo', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--system', type=Path, help='Windows System32 directory when cross-packaging')
    parser.add_argument('--qcustomplot', type=Path, required=True, help='QCustomPlot source directory')
    parser.add_argument('--eigen', type=Path, required=True, help='Eigen source directory')
    parser.add_argument('--qt-source', type=Path, required=True, help='Matching Qt source tar archive (for notices)')
    parser.add_argument('--map-data', type=Path, required=True, help='Packed Natural Earth directory from scripts/pack_natural_earth.py')
    parser.add_argument('--source-ref', required=True, help='Public repository tag or commit matching the executable')
    args = parser.parse_args()
    win = args.platform == 'windows'
    prefix, output = args.prefix.resolve(), args.output.resolve()
    if output.exists():
        raise SystemExit('Output already exists; choose a fresh staging directory.')
    output.mkdir(parents=True)
    owners, records, used = {}, {}, set()
    for path in (prefix / 'conda-meta').glob('*.json'):
        record = json.loads(path.read_text())
        key = path.stem
        records[key] = record
        for file in record.get('files', []):
            owners[file.replace('\\', '/').lower()] = key
    copied = {}
    inventory = []
    queue = deque()
    caches = {}
    for directory, children, _ in os.walk(args.cache):
        directory = Path(directory)
        if (directory / 'info/index.json').exists():
            caches[directory.name] = directory
            children.clear()

    def copy(source, destination):
        source, destination = Path(source), Path(destination)
        identity = str(destination.relative_to(output))
        if identity in copied:
            return
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
        copied[identity] = source
        if source.is_relative_to(prefix):
            relative = source.relative_to(prefix).as_posix().lower()
            owner = owners.get(relative) or owners.get(source.resolve().relative_to(prefix).as_posix().lower())
            if not owner:
                raise RuntimeError(f'No package ownership for {source.name}')
            used.add(owner)
        queue.append((source, destination))

    executable = output / ('EDITools.exe' if win else 'bin/EDITools.bin')
    copy(args.binary.resolve(), executable)
    library = prefix / ('Library/bin' if win else 'lib')
    pluginroot = prefix / ('Library/lib/qt6/plugins' if win else 'lib/qt6/plugins')
    if not pluginroot.exists():
        pluginroot = prefix / ('Library/plugins' if win else 'plugins')
    if win and not (pluginroot / 'styles/qmodernwindowsstyle.dll').exists():
        raise RuntimeError('The Qt Windows 11 style plugin is missing from the build prefix.')
    suffix = '.dll' if win else '.so'
    platforms = ('qwindows.dll', 'qoffscreen.dll') if win else ('libqxcb.so', 'libqoffscreen.so')
    for name in platforms:
        copy(pluginroot / 'platforms' / name, output / 'plugins/platforms' / name)
    for group in ('imageformats', 'iconengines', 'printsupport', 'styles', 'tls') + (() if win else ('xcbglintegrations',)):
        for path in sorted((pluginroot / group).glob('*' + suffix)):
            if group == 'imageformats' and path.stem.removeprefix('lib') not in ('qjpeg', 'qsvg', 'qico', 'qgif'):
                continue
            copy(path, output / 'plugins' / group / path.name)
    candidates = {p.name.lower(): p for folder in (library, prefix) for p in folder.glob('*.dll')} if win else {}
    # Qt loads OpenSSL dynamically; it does not necessarily appear in the import table.
    ssl = ('libssl-3-x64.dll', 'libcrypto-3-x64.dll') if win else ('libssl.so.3', 'libcrypto.so.3')
    for name in ssl:
        path = candidates.get(name.lower()) if win else library / name
        if path and path.exists():
            copy(path, output / (name if win else 'lib/' + name))
    system_dependencies = set()
    while queue:
        source, target = queue.popleft()
        if win:
            import pefile
            pe = pefile.PE(str(source), fast_load=True)
            pe.parse_data_directories(directories=[1, 13])
            names = [entry.dll.decode('ascii') for kind in ('DIRECTORY_ENTRY_IMPORT', 'DIRECTORY_ENTRY_DELAY_IMPORT') for entry in getattr(pe, kind, [])]
            pe.close()
            for name in names:
                dep = candidates.get(name.lower())
                if dep:
                    copy(dep, output / dep.name)
                elif name.lower().startswith(('api-ms-', 'ext-ms-')) or (args.system and (args.system / name).exists()):
                    system_dependencies.add(name)
                else:
                    raise RuntimeError(f'Unresolved DLL {name} required by {source.name}')
        else:
            text = run('ldd', source)
            if 'not found' in text:
                raise RuntimeError(f'Unresolved dependency in {source.name}: {text}')
            paths = re.findall(r'(?:=>\s+|^\s*)(/\S+)\s+\(', text, re.M)
            for name in paths:
                dep = Path(name)
                if re.fullmatch(r'(?:ld-linux.*|lib(?:c|m|dl|pthread|rt|resolv|util)\.so\.\d+)', dep.name):
                    system_dependencies.add(dep.name)
                    continue
                if not dep.is_relative_to(prefix):
                    # Graphics drivers must be supplied by the target Linux installation.
                    if dep.name.startswith(('libGL', 'libEGL', 'libOpenGL')):
                        system_dependencies.add(dep.name)
                        continue
                    raise RuntimeError(f'Unaccounted system dependency: {dep}')
                copy(dep, output / 'lib' / dep.name)
            for needed in run('patchelf', '--print-needed', target).splitlines():
                if '/' in needed:
                    run('patchelf', '--replace-needed', needed, Path(needed).name, target)
            run('patchelf', '--set-rpath', '$ORIGIN/' + os.path.relpath(output / 'lib', target.parent), target)
            run('strip', '--strip-unneeded', target)

    shared = prefix / ('Library/share' if win else 'share')
    for group in ('proj', 'gdal'):
        source = shared / group
        files = [source / 'proj.db', source / 'proj.ini'] if group == 'proj' else [p for p in source.rglob('*') if p.is_file()]
        for path in files:
            if path.exists():
                copy(path, output / 'share' / group / path.relative_to(source))
    certificate = next(p for p in (prefix / 'ssl/cacert.pem', prefix / 'Library/ssl/cacert.pem') if p.exists())
    copy(certificate, output / 'share/certs/cacert.pem')
    if not win:
        for relative in ('etc/fonts', 'share/X11/locale', 'share/X11/xkb'):
            source = prefix / relative
            for path in source.rglob('*'):
                if path.is_file():
                    copy(path, output / relative / path.relative_to(source))
        for path in copied:
            relocate_prefix(output / path, prefix)
    shutil.copy2(args.repo / 'LICENSE', output / 'LICENSE')
    licenses = output / 'licenses'
    licenses.mkdir()
    (licenses / 'Fluent-UI-System-Icons').mkdir()
    for name in ('LICENSE', 'README.md', 'sources.json'):
        shutil.copy2(args.repo / 'resources/icons/fluent' / name, licenses / 'Fluent-UI-System-Icons' / name)
    (licenses / 'QCustomPlot').mkdir()
    shutil.copy2(args.qcustomplot / 'GPL.txt', licenses / 'QCustomPlot/GPL.txt')
    (licenses / 'Eigen').mkdir()
    for path in args.eigen.glob('COPYING*'):
        shutil.copy2(path, licenses / 'Eigen' / path.name)
    for key in sorted(used):
        record = records[key]
        cache = caches.get(key)
        if cache is None:
            raise RuntimeError(f'Package cache is required for license files: {key}')
        files = cache / 'info/licenses'
        if not files.exists() and (cache / 'share/licenses').exists():
            files = cache / 'share/licenses'
        if not files.exists() or not any(files.rglob('*')):
            # Split runtime packages sometimes keep notices in a sibling output.
            import yaml
            recipe = yaml.safe_load((cache / 'info/recipe/meta.yaml').read_text())
            for cached in caches.values():
                sibling = cached / 'info/recipe/meta.yaml'
                if not sibling.exists():
                    continue
                if sibling == cache / 'info/recipe/meta.yaml':
                    continue
                notices = sibling.parent.parent / 'licenses'
                if not notices.exists():
                    continue
                other = yaml.safe_load(sibling.read_text())
                source, candidate = recipe.get('source'), other.get('source')
                same_source = source == candidate or (
                    isinstance(source, dict) and isinstance(candidate, dict)
                    and source.get('url') == candidate.get('url')
                    and recipe['package']['version'] == other['package']['version'])
                if same_source:
                    files = notices
                    break
            else:
                if record.get('license') not in ('Public-Domain', 'blessing'):
                    raise RuntimeError(f'License text missing from package {key}')
                files = None
        if files:
            shutil.copytree(files, licenses / key)
        else:
            notice = 'Public domain. See DEPENDENCIES.json for provenance.\n'
            if record['license'] == 'blessing':
                notice += ('The author disclaims copyright to this source code. In place of\n'
                           'a legal notice, here is a blessing:\n\n'
                           'May you do good and not evil.\n'
                           'May you find forgiveness for yourself and forgive others.\n'
                           'May you share freely, never taking more than you give.\n')
            (licenses / (key + '.txt')).write_text(notice)
        inventory.append({field: record.get(field) for field in ('name', 'version', 'build', 'license', 'url', 'sha256')})
    (output / 'DEPENDENCIES.json').write_text(json.dumps(inventory, indent=2) + '\n')
    with tarfile.open(args.qt_source, 'r|*') as archive:
        for member in archive:
            parts = Path(member.name).parts
            if not member.isfile() or '..' in parts or len(parts) < 3 or parts[1] not in ('qtbase', 'qtsvg'):
                continue
            name = parts[-1].lower()
            notice = name.startswith(('license', 'licence', 'copying', 'copyright', 'readme', 'ftl', 'gpl', 'lgpl'))
            if (len(parts) == 3 and name.startswith('license')) or ('3rdparty' in parts and (notice or name == 'qt_attribution.json')):
                target = licenses / 'Qt-third-party' / Path(*parts[1:])
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(archive.extractfile(member).read())
    copy_map_data(args.map_data, output)
    (output / 'SYSTEM-DEPENDENCIES.txt').write_text('\n'.join(sorted(system_dependencies)) + '\n')
    archive_licenses(output)
    write_launch_note(output, args.platform, args.source_ref)
    (output / ('qt.conf' if win else 'bin/qt.conf')).write_text('[Paths]\nPlugins=' + ('plugins' if win else '../plugins') + '\n')
    if not win:
        launcher = output / 'EDITools'
        launcher.write_text('''#!/bin/sh
set -eu
bundle=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
export LD_LIBRARY_PATH="$bundle/lib"
export QT_PLUGIN_PATH="$bundle/plugins"
export QT_QPA_PLATFORM_PLUGIN_PATH="$bundle/plugins/platforms"
export QT_QPA_PLATFORM="${QT_QPA_PLATFORM:-xcb}"
export FONTCONFIG_FILE="$bundle/etc/fonts/fonts.conf"
export FONTCONFIG_PATH="$bundle/etc/fonts"
export XLOCALEDIR="$bundle/share/X11/locale"
export XKB_CONFIG_ROOT="$bundle/share/X11/xkb"
exec "$bundle/bin/EDITools.bin" "$@"
''')
        launcher.chmod(0o755)
    print(json.dumps({'platform': args.platform, 'files': len(copied), 'packages': len(used), 'output': str(output)}))


if __name__ == '__main__':
    main()
