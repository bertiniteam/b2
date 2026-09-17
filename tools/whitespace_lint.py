#!/usr/bin/env python3
"""Whitespace gate: no tabs, no trailing whitespace, LF line endings, one final newline.

This is the whole of the repository's indentation rule -- spaces, four per level in C++,
Python and CMake -- reduced to what a machine can check without judging style: a tab
character anywhere in a source file is a violation, so are trailing spaces, a CR, and a
file that does not end in exactly one newline.  Indentation DEPTH is not judged; that is a
formatter's job, and this is not a formatter.

Run without arguments to check the repository (the CI gate in .github/workflows/doc_lint.yml
does exactly that; zero tolerance, no baseline).  Run with --fix to rewrite the offending
files in place: tabs expand at 4 columns (str.expandtabs(4), so a tab used for alignment
after the indentation, or in the middle of a line, keeps its column), trailing whitespace
is removed, CRLF becomes LF, and the file ends with one newline.  --fix is the same code
path as the check, so the two can never disagree; a tree that --fix has touched passes.

Paths may be given to check or fix a subset (a file or a directory); the default is the
whole repository from the directory this script lives under.

Exit status: 0 when clean, 1 when violations were found (or, with --fix, when files were
rewritten -- so a pre-commit hook can use it either way).

Covered: .hpp .cpp .h .c .py .pyx .cmake .rst .yml .yaml .toml .sh and every CMakeLists.txt.
Not covered: Markdown (two trailing spaces are a line break there), notebooks, anything
under build/, bld/, _deps/, .git/, python/drafts/.  Pure standard library; no pip deps.
"""

import argparse
import os
import pathlib
import sys

TAB_WIDTH = 4
EXTENSIONS = {'.hpp', '.cpp', '.h', '.c', '.py', '.pyx', '.cmake', '.rst', '.yml', '.yaml', '.toml', '.sh'}
NAMED_FILES = {'CMakeLists.txt'}
SKIP_DIRS = {'build', 'bld', '_deps', '.git', 'drafts', '__pycache__', 'node_modules',
             'bertini_output', 'solve_records', 'cellular_records', 'site'}

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent


def covered(path):
    return path.suffix in EXTENSIONS or path.name in NAMED_FILES


def iter_files(roots):
    for root in roots:
        root = pathlib.Path(root)
        if root.is_file():
            if covered(root):
                yield root
            continue
        for dirpath, dirnames, filenames in os.walk(root):
            dirnames[:] = sorted(d for d in dirnames if d not in SKIP_DIRS)
            for name in sorted(filenames):
                p = pathlib.Path(dirpath) / name
                if covered(p):
                    yield p


def normalize(text):
    """The canonical form of a file's text: LF, no tabs, no trailing whitespace, one final newline."""
    text = text.replace('\r\n', '\n').replace('\r', '\n')
    lines = [line.expandtabs(TAB_WIDTH).rstrip(' \t\f\v') for line in text.split('\n')]
    while lines and lines[-1] == '':
        lines.pop()
    return '\n'.join(lines) + '\n'


def violations(text):
    """Human-readable findings for a file's text, at most a few per kind."""
    found = []
    if '\r' in text:
        found.append('carriage return (CRLF line endings)')
    for number, line in enumerate(text.split('\n'), start=1):
        if '\t' in line:
            found.append(f'line {number}: tab character')
        if line.rstrip('\r') != line.rstrip('\r').rstrip(' \t\f\v'):
            found.append(f'line {number}: trailing whitespace')
        if len(found) >= 6:
            found.append('...')
            break
    if not text.endswith('\n'):
        found.append('no newline at end of file')
    elif text.endswith('\n\n'):
        found.append('more than one newline at end of file')
    return found


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('paths', nargs='*', help='files or directories to check (default: the repository)')
    parser.add_argument('--fix', action='store_true', help='rewrite offending files in place instead of reporting')
    args = parser.parse_args(argv)

    roots = args.paths or [REPO_ROOT]
    bad = 0
    checked = 0
    for path in iter_files(roots):
        checked += 1
        try:
            raw = path.read_bytes()
        except OSError as e:
            print(f'{path}: unreadable ({e})', file=sys.stderr)
            bad += 1
            continue
        if raw == b'':
            continue                                   # an empty file is fine as it is
        text = raw.decode('utf-8', errors='surrogateescape')
        clean = normalize(text)
        if clean == text:
            continue
        bad += 1
        rel = os.path.relpath(path, REPO_ROOT)
        if args.fix:
            path.write_bytes(clean.encode('utf-8', errors='surrogateescape'))
            print(f'fixed  {rel}')
        else:
            print(f'{rel}:')
            for v in violations(text):
                print(f'    {v}')

    verb = 'rewritten' if args.fix else 'with violations'
    print(f'== whitespace lint: {checked} files checked, {bad} {verb} ==')
    if bad and not args.fix:
        print('run  python tools/whitespace_lint.py --fix  to normalize them (see .editorconfig for the rule)')
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
