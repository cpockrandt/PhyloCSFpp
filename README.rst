PhyloCSF++: A fast and user-friendly implementation of PhyloCSF with annotation tools
=====================================================================================

.. image:: https://img.shields.io/conda/dn/bioconda/phylocsfpp.svg?style=flag&label=BioConda%20install
    :target: https://anaconda.org/bioconda/phylocsfpp
    :alt: BioConda Install
.. image:: https://img.shields.io/github/downloads/cpockrandt/phylocsfpp/total.svg
    :target: https://github.com/cpockrandt/PhyloCSFpp/releases/latest
    :alt: Github All Releases
.. image:: https://travis-ci.com/cpockrandt/PhyloCSFpp.svg?branch=master
    :target: https://travis-ci.com/cpockrandt/PhyloCSFpp
    :alt: Travis CI
.. image:: https://img.shields.io/badge/License-AGPLv3-blue.svg
    :target: https://opensource.org/licenses/AGPL-3.0
    :alt: AGPL v3 License

.. contents::
   :local:
   :depth: 2

Introduction
^^^^^^^^^^^^

PhyloCSF can identify protein-coding regions in the genome based on multiple-sequence alignments.
PhyloCSF++ is an implementation of the original methods `PhyloCSF`_ and `PhyloCSF HMM`_.
It allows you to easily create browser tracks for any genome to identify coding regions, score local alignments or
annotate GFF/GTF files with PhyloCSF and its confidence scores.

If you find our implementation useful and use it in your work, please consider citing it:

    Christopher Pockrandt, Martin Steinegger, Steven Salzberg. **PhyloCSF++: A fast and user-friendly implementation of PhyloCSF with annotation tools**. `bioRxiv`_, 2021.

Please also consider citing the original method papers:

    Lin MF et al. PhyloCSF: a comparative genomics method to distinguish protein-coding and non-coding regions. Bioinformatics, 2011.

    Mudge JM et al. Discovery of high-confidence human protein-coding genes and exons by whole-genome PhyloCSF helps elucidate 118 GWAS loci. Genome Research, 2019.

.. _bioRxiv: https://doi.org/10.1101/2021.03.10.434297

Installation
^^^^^^^^^^^^

Binaries
""""""""

Download the `latest release <https://github.com/cpockrandt/PhyloCSFpp/releases/latest>`_ (static binary for Linux 64bit).

Bioconda
""""""""

::

    $ conda install -c conda-forge -c bioconda phylocsfpp

Building from source
""""""""""""""""""""

If you want to build it from source, we recommend cloning the git repository as shown below.
You will need the GNU scientific library (gsl), which should be available in your package manager.

::

    $ git clone https://github.com/cpockrandt/PhyloCSFpp.git PhyloCSF++
    $ mkdir PhyloCSF++/build && cd PhyloCSF++/build
    $ cmake .. -DCMAKE_BUILD_TYPE=Release
    $ make

You can then install and run PhyloCSF++ as follows

::

    $ sudo make install
    $ phylocsf++ --help

or run the binary directly:

::

    $ ./phylocsf++ --help

**Requirements**

Operating System
  GNU/Linux, Mac

Compiler
  GCC ≥ 4.9, Clang ≥ 3.8

Build system
  CMake ≥ 3.2

Dependencies
  GNU scientific library (gsl), OpenMP, zlib

Getting started
^^^^^^^^^^^^^^^

A detailed list of arguments and explanations can be retrieved with ``--help``:

::

    $ phylocsf++ --help
    $ phylocsf++ build-tracks --help

Please also check out the `FAQ in our wiki <https://github.com/cpockrandt/PhyloCSFpp/wiki>`_.

Building tracks
"""""""""""""""

To produce PhyloCSF tracks including the power track (confidence scores), you can run the following command.
Afterwards, you can use `wigToBigWig <http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/>`_ to transform the outputted ``wig`` files to ``bw`` files and load them into a genome browser.

::

    $ phylocsf++ build-tracks 58mammals hg38.100way.maf

We recommend to additionally smoothen the tracks (i.e., compute posterior probabilities) by setting ``--output-phylo 1``.
This requires to pass the total genome length and a list of known coding regions as training data.
Here we extract the known coding regions from a GFF file using awk (for more information on the file format, please have a look at the `wiki <https://github.com/cpockrandt/PhyloCSFpp/wiki>`_.

::

    $ awk -F'\t' 'BEGIN { OFS="\t" } ($3 == "CDS") { print $1, $7, $8, $4, $5 }' gene_catalogue.gff > CodingExons.txt
    $ phylocsf++ build-tracks --output-phylo 1 --genome-length 3272116950 --coding-exons CodingExons.txt 58mammals hg38.100way.maf

If not all of the species from the model are present in your alignments, reducing the model can speed up the computation significantly. For this pass all species from your alignment, e.g., ``--species Human,Chimp,Gorilla``.

Here is a minimal example with all files at hand in the repository.
It includes a small set of alignments for chicken (galGal6).
The wig files are written to ``./galgal6-tracks``.

::

    $ gunzip example/galGal6_chr22_25_28_each_30k_bases.maf.gz
    $ phylocsf++ build-tracks \
        --threads 4 --output ./galgal6-tracks \
        --output-phylo 1 --genome-length 1065365434 --coding-exons example/galGal6_coding_exons.txt \
        53birds example/galGal6_chr22_25_28_each_30k_bases.maf

Annotating GFF files with tracks
""""""""""""""""""""""""""""""""

If you have tracks in bigWig computed or downloaded, PhyloCSF++ can annotate CDS with PhyloCSF and confidence scores:

::

    $ phylocsf++ annotate-with-tracks /path/to/PhyloCSF+1.bw genes.gff

For this you need to have all six files in the same directory (PhyloCSF+1.bw, PhyloCSF+2.bw, etc.) as well as PhyloCSFpower.bw if you also want to compute confidence scores.

Here is a minimal example with all files at hand in the repository.
It includes a few transcripts from chicken (galGal6) and precomputed tracks.
The output is written into the same directory with the suffix ``.PhyloCSF++.gtf``.

::

    $ phylocsf++ annotate-with-tracks ./example/tracks/PhyloCSF+1.bw ./example/galGal6_chr22_25_28_subset_ensGene.gtf
    $ less ./example/galGal6_chr22_25_28_subset_ensGene.PhyloCSF++.gtf

For some species you can download complete track files either on the `Broad institute server <https://data.broadinstitute.org/compbio1/PhyloCSFtracks/>`_ or here: ftp://ftp.ccb.jhu.edu/pub/software/phylocsfpp/

Scoring alignments
""""""""""""""""""

If you want to score alignments and not create tracks for an entire genome, simply run:

::

    $ phylocsf++ score-msa 58mammals hg38.100way.maf

You can also specify the strategy (fixed, mle and omega; default: mle) and choose which scores to compute (PhyloCSF score, ancestral sequence composition score, branch length score).

NOTE: compared to the original implementation of PhyloCSF, PhyloCSF++ only scores the forward strand starting from the first base.
For other frames and strands, you need to remove the first 1-2 bases and/or compute the reverse complement of the sequences.

To make these scores easier to interpret, we added the mode ``fixed_mean``.
It scores every codon in the MSA, computes posterior probabilities and computes a mean over all codons.
Hence, the final score is in the interval [-15, +15] just as the tracks.

Annotating GFF files with MMseqs
""""""""""""""""""""""""""""""""

If you don't have tracks available for your genome of interest, PhyloCSF++ can annotate CDS with PhyloCSF and confidence scores by computing an alignment on the fly using MMseqs:

::

    $ phylocsf++ annotate-with-mmseqs genomes.txt 58mammals genes.gff

``genomes.txt`` has to contain the paths to genomes from the selected model to align to (one per line).
After all CDS lines were extracted and aligned with MMseqs, PhyloCSF++ scores each CDS alignment with the sub-tool ``score-msa``.

Motivation
^^^^^^^^^^

We think that PhyloCSF is a very useful method for gene finding and annotation.
Unfortunately no binaries are available and we think the outdated Ocaml code might be difficult to get running for inexperienced users.
To build tracks the user also has to set up their own pipeline and do some coding.
Hence, we thought it would be helpful to make an easy-to-use program that merges all necessary steps into a single step to quickly create tracks for entire genomes.
As part of this project we computed tracks for more species and included them into the UCSC genome browser as well as offer them for download:

ftp://ftp.ccb.jhu.edu/pub/software/phylocsfpp

License
^^^^^^^

This is an implementation of the original methods (`PhyloCSF`_ and `PhyloCSF HMM`_), which were released under the GNU AGPL v3 and Apache License v2.
We have reimplemented the core algorithms (originally written in OCaml and Python) in C++, they were not changed except for running time improvements or where explicitly stated in the source code.

.. _PhyloCSF: https://github.com/mlin/PhyloCSF
.. _PhyloCSF HMM: https://github.com/iljungr/PhyloCSFCandidateCodingRegions