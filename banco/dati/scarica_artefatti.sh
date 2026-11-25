#!/bin/bash
# Scarica e scompatta gli artefatti dell'ultimo giro di un flusso di giu-cloud/lens.
#   scarica_artefatti.sh <flusso.yml> <cartella> [id del giro]
flusso=$1; out=$2; run=$3
mkdir -p "$out"
[ -z "$run" ] && run=$(gh api "repos/giu-cloud/lens/actions/workflows/$flusso/runs?per_page=1" --jq '.workflow_runs[0].id')
gh api "repos/giu-cloud/lens/actions/runs/$run/artifacts?per_page=100" --jq '.artifacts[] | "\(.id) \(.name)"' | while read id n; do
  gh api "repos/giu-cloud/lens/actions/artifacts/$id/zip" > /tmp/art_$id.zip && python3 -c "import zipfile,sys; zipfile.ZipFile(sys.argv[1]).extractall(sys.argv[2])" /tmp/art_$id.zip "$out" && rm -f /tmp/art_$id.zip
  echo "$n scompattato"
done
