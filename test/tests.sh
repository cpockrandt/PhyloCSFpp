#!/bin/bash -ex

PhyloCSFpp="build/phylocsf++"
THREADS=$1

### build-tracks

# download data: wget -m -nd ftp://ftp.ccb.jhu.edu/pub/pocki/phylocsf++/tests/s
gunzip test/chr22_25_28_each_30k_bases.maf.gz

$PhyloCSFpp build-tracks --genome-length 1065365434 --coding-exons test/ChickenCodingExons.txt \
  --output-phylo 1 --output-raw-phylo 1 --output-power 1 --threads $THREADS --output ./build/results \
  53birds test/chr22_25_28_each_30k_bases.maf

diff -r ./build/results ./test/expected_results

### score-msa

species=""
# --species Lizard,Dolphin,Orangutan,Platypus,Opossum,Weddell_seal,Human,Gibbon,Stickleback,Baboon,Chimp
# --species anoCar,gasAcu,hg,lepWed,monDom,nomLeu,ornAna,panTro,papHam,ponAbe,turTru

echo "FIXED mode"
$PhyloCSFpp score-msa --threads $THREADS --strategy fixed --comp-phylo 1 --comp-anc 1 --comp-bls 1 $species 100vertebrates example/small/chr22.50alignments.maf # no randomization
DIFF_FIXED=$(diff <(awk '!($0 ~ /^#/)' example/small/chr22.50alignments.maf.scores) example/small/PhyloCSFpp-results/chr22.50alignments.fixed.scores) # remove header line
echo -e "${DIFF_FIXED}"

echo "MLE mode"
$PhyloCSFpp score-msa --threads $THREADS --strategy mle --comp-phylo 1 --comp-anc 1 --comp-bls 1 100vertebrates example/small/chr22.50alignments.maf
DIFF_MLE=$(awk -F$'\t' '($1 != $8 || $2 != $9 || $3 != $10 || $4 != $11 || ($5-$12)^2 > 0.001 || ($6-$13)^2 > 0.001 || $7 != $14)' <(paste <(awk '!($0 ~ /^#/)' example/small/chr22.50alignments.maf.scores) example/small/PhyloCSFpp-results/chr22.50alignments.mle.scores))
echo -e "${DIFF_MLE}"

echo "OMEGA mode"
$PhyloCSFpp score-msa --threads $THREADS --strategy omega --comp-phylo 1 --comp-anc 0 --comp-bls 1 100vertebrates example/small/chr22.50alignments.maf
DIFF_OMEGA=$(awk -F$'\t' '($1 != $7 || $2 != $8 || $3 != $9 || $4 != $10 || ($5-$11)^2 > 0.1 || $6 != $12)' <(paste <(awk '!($0 ~ /^#/)' example/small/chr22.50alignments.maf.scores) example/small/PhyloCSFpp-results/chr22.50alignments.omega.scores))
#DIFF_OMEGA=$(diff example/small/chr22.50alignments.maf.scores example/small/PhyloCSFpp-results/chr22.50alignments.omega.scores)
echo -e "${DIFF_OMEGA}"

if [ -n "${DIFF_FIXED}" ] || [ -n "${DIFF_MLE}" ] || [ -n "${DIFF_OMEGA}" ]; then # is any diff non-empty?
  exit 1
fi