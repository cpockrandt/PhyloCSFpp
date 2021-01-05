#!/bin/bash -e

MAF="/home/chris/Downloads/chr22.head.maf"
FASTA="/home/chris/Downloads/chr22.head.fasta"
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
