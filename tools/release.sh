#!/usr/bin/env bash
# Make a GitHub release, but only after checking it matches the program it ships.
# Zenodo archives every release the moment it is published, so the check has to
# happen here, before publishing, not afterwards.
#
#   tools/release.sh v3.1.0 "What changed, in a sentence or two"
set -euo pipefail
tag=${1:?usage: tools/release.sh vX.Y.Z "release notes"}
notes=${2:?usage: tools/release.sh vX.Y.Z "release notes"}
cd "$(dirname "$0")/.."

if [ -n "$(git status --porcelain)" ]; then
  echo "REFUSED: there are uncommitted changes. Commit them first."; exit 1
fi
git fetch -q --tags origin
if [ "$(git rev-parse HEAD)" != "$(git rev-parse '@{u}')" ]; then
  echo "REFUSED: this copy and GitHub differ. Push (or pull) first."; exit 1
fi
if git rev-parse -q --verify "refs/tags/$tag" >/dev/null; then
  echo "REFUSED: tag $tag already exists."; exit 1
fi

python3 tools/check_version.py --tag "$tag"

gh release create "$tag" --title "RaGCAn $tag" --notes "$notes" --target "$(git rev-parse HEAD)"
echo "Released $tag. Zenodo will mint its DOI within a minute or two."
