PhyloCSF++: A fast and easy-to-use implementation of PhyloCSF
=============================================================

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
It allows you to easily create browser tracks for any genome to identify coding regions or score local alignments.

If you find our implementation useful, please consider citing it:

Christopher Pockrandt, Martin Steinegger, Steven Salzberg. **PhyloCSF++: Making PhyloCSF available for more species to improve gene annotation pipelines**. `Bioinformatics`_, 2021.

Please also consider citing the original method paper:

Lin MF, Jungreis I, and Kellis M. PhyloCSF: a comparative genomics method to distinguish protein-coding and non-coding regions. Bioinformatics, 2011.

.. _Bioinformatics: https://doi.org/10.1093/bioinformatics/btaa222

Installation
^^^^^^^^^^^^

Binaries
""""""""

Download the `latest release <https://github.com/cpockrandt/PhyloCSFpp/releases/latest>`_ (static binary for Linux 64bit).

Bioconda
""""""""

::

    $ conda install -c bioconda phylocsfpp

Building from source
""""""""""""""""""""

If you want to build it from source, we recommend cloning the git repository as shown below.
Please note that building from source can easily take 10 minutes and longer depending on your machine and compiler.

::

    $ git clone https://github.com/cpockrandt/PhyloCSFpp.git PhyloCSF++
    $ mkdir PhyloCSF++/build && cd PhyloCSF++/build
    $ cmake .. -DCMAKE_BUILD_TYPE=Release
    $ make

You can then install and run PhyloCSF++ as follows

::

    $ sudo make install
    $ phylocsf++

or run the binary directly:

::

    $ ./phylocsf++

**Requirements**

Operating System
  GNU/Linux, Mac

Compiler
  GCC ≥ 4.9, LLVM/Clang ≥ 3.8 (requires C++11)

Build system
  CMake ≥ 2.8.12

Dependencies
  OpenMP, GNU scientific library (gsl)

Getting started
^^^^^^^^^^^^^^^

A detailed list of arguments and explanations can be retrieved with ``--help``:

::

    $ phylocsf++ --help
    $ phylocsf++ tracks --help
    $ phylocsf++ regions --help

Please also check out the `FAQ in our wiki <https://github.com/cpockrandt/PhyloCSFpp/wiki>`_.

Building tracks
"""""""""""""""

To produce PhyloCSF tracks including the Power track (confidence scores), you can run the following command.
Afterwards, you can use `wigToBigWig <http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/>`_ to transform the outputted ``wig`` files to ``bw`` files and load them into a genome browser.

::

    $ phylocsf++ tracks 58mammals hg38.100way.maf

To additionally smoothen the tracks, you can set ``--output-phylo 1``.
This requires to pass the total genome length and a list of known coding regions (check the `wiki <https://github.com/cpockrandt/PhyloCSFpp/wiki>`_ for a command how to extract it from a gtf/gff file) as training data.

::

    $ phylocsf++ tracks --output-phylo 1 --genome-length 3272116950 --coding-exons HumanCodingExonsV29.txt 58mammals hg38.100way.maf

Scoring alignments
""""""""""""""""""

If you want to score entire alignments and not every codon separately, simply run:

::

    $ phylocsf++ region 100vertebrates hg38.100way.maf

You can also specify the strategy (default: mle) and choose which scores to compute (PhyloCSF, ancestral sequence composition score, branch length score).

Training my own model
"""""""""""""""""""""

For training your own model, a phylogenetic tree with evolutionary distances, as well as codon frequencies and codon substitution rates for both coding and non-coding regions are required.
At the moment neither PhyloCSF, nor PhyloCSF++ has a tool to compute this model, but we are planning to include it into PhyloCSF++ in the near future.

Motivation
^^^^^^^^^^

We think that PhyloCSF is a very useful method for gene finding and annotation.
Unfortunately no binaries are available and we think the outdated Ocaml code might be difficult to get running for inexperienced users.
To build tracks the user also has to set up their own pipeline and do some coding.
Hence, we thought it would be helpful to make an easy-to-use program that merges all necessary steps into a single step to quickly create tracks for entire genomes.
As part of this project we computed tracks for many algorithms and included them into the UCSC genome browser as well as offer them for download:

ftp://ftp.ccb.jhu.edu/pub/pocki/phylocsf++

License
^^^^^^^

This is an implementation of the original methods (`PhyloCSF`_ and `PhyloCSF HMM`_), which were released under the GNU AGPL v3 and Apache License v2.
We have reimplemented the core algorithms (originally written in OCaml and Python) in C++, they were not changed except for running time improvements or where explicitly stated in the source code.

.. _PhyloCSF: https://github.com/mlin/PhyloCSF
.. _PhyloCSF HMM: https://github.com/iljungr/PhyloCSFCandidateCodingRegions