#!/bin/bash -e

MAF="/home/chris/Downloads/chr22.maf"
FASTA="/home/chris/Downloads/chr22.fasta"
NAMES="/home/chris/dev-uni/PhyloCSF++/commonNames_assemblies.txt"
PHYLOCSF="/home/chris/dev-uni/PhyloCSF/PhyloCSF"

# 1. MAF to FASTA (and changed names)
# awk -F$'\t' 'FNR==NR { f[$2]=$1; next } { if ($0 ~ /^>/) { $0 = ">" f[substr($0, 2)]; } print }' $NAMES \
# 	<(awk -F' ' '($0 ~ /^a/) { print "" } ($0 ~ /^s/) { print ">" substr($2, 1, index($2, ".")-1) "\n" $7 }' $MAF | tail -n+2) > $FASTA

# awk '($0 == "") { id++; next } (id == 1-)'
I=0
while true; do
	# echo "$I"

	awk -v ID="$I" '($0 == "") { i++; next } (i == ID) { print } (i > ID) { exit }' $FASTA > /tmp/phylo.fasta

	echo -en "$I\tFixed:\t"
	$PHYLOCSF --strategy=fixed --frames=1 --bls --ancComp --removeRefGaps 100vertebrates /tmp/phylo.fasta

	echo -en "$I\tMLE:\t"
	$PHYLOCSF --strategy=mle --frames=1 --bls --ancComp --removeRefGaps 100vertebrates /tmp/phylo.fasta

	echo -en "$I\tOmega:\t"
	$PHYLOCSF --strategy=omega --frames=1 --bls --ancComp --removeRefGaps 100vertebrates /tmp/phylo.fasta

	echo ""

	I=$((I+1))
done

# format this file's output with:
# awk 'BEGIN{ OFS="\t"; print "ALN-ID", "FIXED", "FIXED-ANC", "MLE", "MLE-ANC", "OMEGA", "BLS" } ($0 != "" && !($0 ~ /^\/tmp/)) { if ($3 ~ /^\/tmp/) { $3 = "nan" } if ($2 == "Fixed:") { bls[$1] = $4 } anc[$1][$2] = $5; f[$1][$2] = $3 } END { for (id in f) print id, f[id]["Fixed:"], anc[id]["Fixed:"], f[id]["MLE:"], anc[id]["MLE:"], f[id]["Omega:"], bls[id] }' chr22.results > chr22.results.formatted
