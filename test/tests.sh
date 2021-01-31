#!/bin/bash -ex

PhyloCSFpp="build/phylocsf++"

$PhyloCSFpp region --threads 6 --strategy fixed --comp-phylo 1 --comp-anc 1 --comp-bls 1 100vertebrates example/chr22.516alignments.maf
DIFF_FIXED=$(awk -F$'\t' '($1 != $8 || $2 != $9 || $3 != $10 || $4 != $11 || ($5-$12)^2 > 0.00000001 || ($6-$13)^2 > 0.00000001 || ($7-$14)^2 > 0.00000001)' <(paste example/chr22.516alignments.maf.scores example/chr22.516alignments.maf.PhyloCSF.fixed.scores))
echo -e "${DIFF_FIXED}"

$PhyloCSFpp region --threads 6 --strategy mle --comp-phylo 1 --comp-anc 1 --comp-bls 1 100vertebrates example/chr22.516alignments.maf
DIFF_MLE=$(awk -F$'\t' '($1 != $8 || $2 != $9 || $3 != $10 || $4 != $11 || ($5-$12)^2 > 0.00000001 || ($6-$13)^2 > 0.00000001 || ($7-$14)^2 > 0.00000001)' <(paste example/chr22.516alignments.maf.scores example/chr22.516alignments.maf.PhyloCSF.mle.scores))
echo -e "${DIFF_MLE}"

$PhyloCSFpp region --threads 6 --strategy omega --comp-phylo 1 --comp-anc 0 --comp-bls 1 100vertebrates example/chr22.516alignments.maf
DIFF_OMEGA=$(awk -F$'\t' '($1 != $7 || $2 != $8 || $3 != $9 || $4 != $10 || ($5-$11)^2 > 0.1 || ($6-$12)^2 > 0.00000001)' <(paste example/chr22.516alignments.maf.PhyloCSF.omega.scores example/chr22.516alignments.maf.omega.scores))
echo -e "${DIFF_OMEGA}"

if [ -n "${DIFF_FIXED}" ] || [ -n "${DIFF_MLE}" ] || [ -n "${DIFF_OMEGA}" ]; then # is any diff non-empty?
  exit 1
fi