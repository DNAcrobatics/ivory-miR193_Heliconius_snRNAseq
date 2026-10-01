#!/usr/bin/env python3
"""New intake utility: collect code snapshots and inventory other project files.

Python 3 standard library only. Does not run analyses or change source files.
"""
import argparse
import csv
import fnmatch
import hashlib
import io
import json
import os
from pathlib import Path
import tarfile
import tempfile

SUFFIXES = {'.r', '.py', '.sh', '.bash', '.slurm', '.sbatch', '.rmd', '.qmd',
            '.ipynb', '.jl', '.pl', '.awk', '.nf', '.smk', '.snakefile',
            '.yaml', '.yml', '.toml', '.ini', '.cfg', '.json', '.wdl'}
NAMES = {'makefile', 'snakefile', 'dockerfile', 'containerfile',
         'requirements.txt', 'renv.lock', 'description', 'namespace',
         'sessioninfo.txt', 'session_info.txt', 'software_versions.txt',
         'versions.txt', 'module_list.txt', 'pip_freeze.txt', 'conda_explicit.txt'}
SKIP_DIRS = {'.git', '.svn', '.venv', 'venv', 'node_modules', '__pycache__',
             '.cache', '.ipynb_checkpoints', '.Rproj.user', 'site-packages', 'pkgs'}
EXCLUDE = {'.env', 'credentials.json', 'service-account.json', 'id_rsa', 'id_ed25519'}
FIELDS = ['root_label', 'original_path', 'relative_path', 'size_bytes',
          'mtime_ns', 'included', 'sha256', 'reason', 'archive_path']


def add_bytes(archive, name, payload):
    info = tarfile.TarInfo(name)
    info.size = len(payload)
    info.mode = 0o644
    archive.addfile(info, io.BytesIO(payload))


def collect(roots, output, max_bytes, globs):
    output = Path(output).absolute()
    if output.exists():
        raise FileExistsError('Choose a new --out path: ' + str(output))
    output.parent.mkdir(parents=True, exist_ok=True)
    inventory, statuses = [], []
    fd, tmp_name = tempfile.mkstemp(prefix='.paper_code_', dir=output.parent)
    os.close(fd)
    temporary = Path(tmp_name)
    try:
        with tarfile.open(temporary, 'w:gz') as archive:
            for index, supplied in enumerate(roots, 1):
                root = Path(supplied).expanduser().absolute()
                label = '{:02d}_{}'.format(index, root.name or 'root')
                status = {'root_label': label, 'path': str(root),
                          'status': 'scanned', 'walk_errors': []}
                statuses.append(status)
                if not root.is_dir() or root.is_symlink():
                    status['status'] = 'missing_or_not_directory_or_symlink'
                    continue
                def onerror(error):
                    status['walk_errors'].append(str(error))
                for folder, dirs, names in os.walk(root, followlinks=False, onerror=onerror):
                    dirs[:] = sorted(d for d in dirs if d not in SKIP_DIRS
                                     and not (Path(folder) / d).is_symlink())
                    # Preserve renv.lock/config but omit installed renv libraries.
                    if Path(folder).name == 'renv':
                        dirs[:] = [d for d in dirs if d not in {'library', 'staging', 'cache'}]
                    for name in sorted(names):
                        path = Path(folder) / name
                        if path.absolute() in {output, temporary}:
                            continue
                        relative = path.relative_to(root).as_posix()
                        row = dict.fromkeys(FIELDS, '')
                        row.update(root_label=label, original_path=str(path),
                                   relative_path=relative, included='no')
                        inventory.append(row)
                        try:
                            if path.is_symlink():
                                row['reason'] = 'symlink_not_followed'
                                continue
                            stat = path.stat()
                            row.update(size_bytes=stat.st_size, mtime_ns=stat.st_mtime_ns)
                            if not path.is_file():
                                row['reason'] = 'not_regular_file'
                                continue
                            if name.lower() in EXCLUDE:
                                row['reason'] = 'excluded_filename'
                                continue
                            selected = (path.suffix.lower() in SUFFIXES or name.lower() in NAMES
                                        or any(fnmatch.fnmatch(relative, g) for g in globs))
                            if not selected:
                                row['reason'] = 'metadata_only_noncode'
                                continue
                            if stat.st_size > max_bytes:
                                row['reason'] = 'selected_but_over_size_cap'
                                continue
                            with path.open('rb') as handle:
                                payload = handle.read(max_bytes + 1)
                            after = path.stat()
                            if len(payload) > max_bytes:
                                row['reason'] = 'selected_but_grew_over_size_cap'
                                continue
                            if (stat.st_size, stat.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
                                row['reason'] = 'changed_during_read_retry_required'
                                continue
                            archived = 'source_collection/{}/{}'.format(label, relative)
                            add_bytes(archive, archived, payload)
                            row.update(included='yes', sha256=hashlib.sha256(payload).hexdigest(),
                                       reason='included_code_or_explicit_glob', archive_path=archived)
                        except OSError as error:
                            row['reason'] = 'read_error: ' + str(error)
            buffer = io.StringIO(newline='')
            writer = csv.DictWriter(buffer, FIELDS, delimiter='\t')
            writer.writeheader()
            writer.writerows(inventory)
            add_bytes(archive, 'source_collection/source_inventory.tsv', buffer.getvalue().encode())
            summary = {'roots': statuses, 'inventoried_files': len(inventory),
                       'included_files': sum(r['included'] == 'yes' for r in inventory),
                       'max_file_bytes': max_bytes, 'include_globs': globs,
                       'skipped_directory_names': sorted(SKIP_DIRS),
                       'limitations': ['Symlinks and excluded directory trees are not followed.',
                                       'Noncode content is omitted unless explicitly selected.',
                                       'This is intake material, not a verified paper release.']}
            add_bytes(archive, 'source_collection/collection_summary.json',
                      (json.dumps(summary, indent=2) + '\n').encode())
            add_bytes(archive, 'source_collection/README.txt',
                      b'Original code snapshots for ivory:mir-193 paper repository intake.\n'
                      b'Consult inventory and summary for provenance and missing/skipped files.\n'
                      b'No analyses were run.\n')
        # Atomic publication with no overwrite, including concurrent invocations.
        os.link(temporary, output)
    finally:
        temporary.unlink(missing_ok=True)
    print(json.dumps({'output': str(output), 'included_files': summary['included_files'],
                      'inventoried_files': len(inventory), 'roots': statuses}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', action='append', required=True, help='Repeat for each project')
    parser.add_argument('--out', required=True, help='New .tar.gz path; refuses overwrite')
    parser.add_argument('--max-file-mib', type=int, default=20, help='Default 20 MiB per file')
    parser.add_argument('--include-glob', action='append', default=[], help='Additional relative-path glob')
    args = parser.parse_args()
    if args.max_file_mib < 1:
        raise ValueError('--max-file-mib must be at least 1')
    collect(args.root, args.out, args.max_file_mib * 1024 * 1024, args.include_glob)


if __name__ == '__main__':
    main()
