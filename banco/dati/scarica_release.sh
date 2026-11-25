#!/bin/bash
# Scarica da giu-cloud/lens le release che cominciano con un prefisso: l'archivio .tar se c'è, altrimenti i file uno per uno.
#   scarica_release.sh <prefisso del tag> <cartella> [rinomina_da rinomina_a]
pre=$1; out=$2; da=$3; a=$4
mkdir -p "$out"
for tag in $(gh api "repos/giu-cloud/lens/releases?per_page=100" --jq ".[] | select(.tag_name|startswith(\"$pre\")) | .tag_name"); do
  gh api "repos/giu-cloud/lens/releases/tags/$tag" --jq '.assets[] | "\(.id) \(.name)"' > /tmp/asset_$tag.txt
  t=$(grep "\.tar$" /tmp/asset_$tag.txt | head -1)
  if [ -n "$t" ]; then
    set -- $t; gh api -H "Accept: application/octet-stream" repos/giu-cloud/lens/releases/assets/$1 > /tmp/$2 && tar xf /tmp/$2 -C "$out" && rm -f /tmp/$2
  else
    grep "\.grb2$" /tmp/asset_$tag.txt | while read id n; do
      d=$n; [ -n "$da" ] && d=${n/$da/$a}
      [ -s "$out/$d" ] || gh api -H "Accept: application/octet-stream" repos/giu-cloud/lens/releases/assets/$id > "$out/$d"
    done
  fi
  echo "$tag fatto"
done
