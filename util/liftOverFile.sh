#!/bin/bash -ex

# Source: https://www.biostars.org/p/81185

IN_FOLDER=$1
OUT_FOLDER=$2
CHAIN=$3
CHROM_SIZES=$4

for d in $(ls ${IN_FOLDER}/*.bw); do
  f=$(basename $d)

  IN=$d
  OUT="${OUT_FOLDER}/${f}"

  if [ -f "${OUT}" ]; then
    continue
  fi

  IN_BEDGRAPH="${IN}.bedgraph"
  OUT_BEDGRAPH="${OUT}.bedgraph"
  OUT_BEDGRAPH_SORTED_UNIQUE="${OUT}.sorted.unique.bedgraph"
  OUT_UNMAPPED="${OUT}.unmapped"

  bigWigToBedGraph $IN $IN_BEDGRAPH
  liftOver $IN_BEDGRAPH $CHAIN $OUT_BEDGRAPH $OUT_UNMAPPED
  bedtools sort -i $OUT_BEDGRAPH | awk -F$'\t' 'BEGIN{OFS="\t"} { if (NR>1) { if ($2+0<p3+0 && p1 == $1) { p3=$2 } print p1, p2, p3, p4 } p1=$1; p2=$2; p3=$3; p4=$4; } END{ print p1, p2, p3, p4 }' > $OUT_BEDGRAPH_SORTED_UNIQUE
  bedGraphToBigWig $OUT_BEDGRAPH_SORTED_UNIQUE $CHROM_SIZES $OUT

  awk -F$'\t' '!($0 ~ /^#/) && $4 != "0" { s+=$3-$2 } END { printf "%d, %d\n", s, s/3}' $OUT_UNMAPPED

  rm $IN_BEDGRAPH $OUT_BEDGRAPH $OUT_UNMAPPED $OUT_BEDGRAPH_SORTED_UNIQUE
done
