#!/usr/bin/env python3
"""Check that a release matches the program it ships.

Always checked:
    the version in RaGCAn.py (__version__) equals the version in CITATION.cff.

With --tag vX.Y.Z, also checked:
    the tag has the form vX.Y.Z and is newer than every earlier release tag;
    if RaGCAn.py has changed since the previous release, the tag equals the version
    inside RaGCAn.py. A release that changes only supporting files may keep the old
    program version, because the program itself is the same file.

Exit status 0 when everything agrees, 1 otherwise. Used by tools/release.sh before a
release is created, and by the GitHub check afterwards.
"""
import argparse
import pathlib
import re
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
TAG = re.compile(r'v(\d+)\.(\d+)\.(\d+)')


def program_version():
    m = re.search(r"^__version__\s*=\s*['\"]([^'\"]+)['\"]", (ROOT / 'RaGCAn.py').read_text(), re.M)
    return m.group(1) if m else None


def citation_version():
    m = re.search(r'^version:\s*"?([^"\s]+)"?\s*$', (ROOT / 'CITATION.cff').read_text(), re.M)
    return m.group(1) if m else None


def git(*args):
    return subprocess.run(['git', *args], cwd=ROOT, capture_output=True, text=True)


def as_tuple(tag):
    return tuple(int(x) for x in TAG.fullmatch(tag).groups())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--tag', help='the release tag about to be created, e.g. v3.1.0')
    tag = ap.parse_args().tag

    errors, notes = [], []
    pv, cv = program_version(), citation_version()
    print(f'RaGCAn.py reports   {pv}')
    print(f'CITATION.cff says   {cv}')
    if pv is None:
        errors.append('No __version__ found in RaGCAn.py.')
    if pv != cv:
        errors.append(f'RaGCAn.py says {pv} but CITATION.cff says {cv}. Make them the same.')

    if tag:
        print(f'release tag         {tag}')
        if not TAG.fullmatch(tag):
            errors.append(f'Tag {tag} is not of the form vX.Y.Z.')
        else:
            earlier = [t for t in git('tag', '--list', 'v*').stdout.split()
                       if TAG.fullmatch(t) and t != tag]
            prev = max(earlier, key=as_tuple) if earlier else None
            print(f'previous release    {prev or "(none found; run git fetch --tags if this is wrong)"}')
            if prev and as_tuple(tag) <= as_tuple(prev):
                errors.append(f'Tag {tag} is not newer than the previous release {prev}.')
            changed = prev is None or git('diff', '--quiet', prev, 'HEAD', '--', 'RaGCAn.py').returncode != 0
            if changed and tag[1:] != pv:
                errors.append(
                    f'RaGCAn.py has changed since {prev or "the start"} but still reports {pv}. '
                    f'Set __version__ in RaGCAn.py and version in CITATION.cff to {tag[1:]} '
                    f'before releasing {tag}.')
            elif not changed and tag[1:] != pv:
                notes.append(f'RaGCAn.py is unchanged since {prev}, so it may keep reporting {pv}. '
                             'This release changes supporting files only.')

    for n in notes:
        print('note: ' + n)
    for e in errors:
        print('ERROR: ' + e)
    print('OK' if not errors else 'REFUSED')
    return 1 if errors else 0


if __name__ == '__main__':
    sys.exit(main())
