#!/usr/bin/env bash
set -euo pipefail
ROOT=/home/asset/projects/IBS_A_OFFICIAL_RERUN_20261105
cd "$ROOT"
mkdir -p data/reference/gwas/archives
FETCH=/home/asset/miniconda3/bin/curl
fetch() {
 url="$1";target="$2"
 if test -f "$target"; then python3 -c 'import sys,zipfile;z=zipfile.ZipFile(sys.argv[1]);assert z.testzip() is None;print("VALID_EXISTING_OFFICIAL_DOWNLOAD",sys.argv[1])' "$target";return;fi
 "$FETCH" --fail --location --connect-timeout 25 --max-time 3600 --retry 2 --retry-delay 5 --dump-header "$target.http_headers.txt" --output "$target.part" "$url"
 python3 -c 'import sys,zipfile;z=zipfile.ZipFile(sys.argv[1]);assert z.testzip() is None;print(z.namelist()[:15])' "$target.part"
 mv -- "$target.part" "$target"
}
case "${1:?static or eur}" in
 static) fetch https://vu.data.surf.nl/index.php/s/lxDgt2dNdNr6DYt/download data/reference/gwas/archives/magma_v1.10_static_official.zip ;;
 eur) fetch https://vu.data.surf.nl/index.php/s/VZNByNwpD8qqINe/download data/reference/gwas/archives/g1000_eur_official.zip ;;
 *) exit 2 ;;
esac
