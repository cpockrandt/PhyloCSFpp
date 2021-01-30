PhyloCSF++: A fast and easy-to-use implementation of PhyloCSF
=============================================================

.. image:: https://img.shields.io/conda/dn/bioconda/phylocsfpp.svg?style=flag&label=BioConda%20install
    :target: https://anaconda.org/bioconda/phylocsfpp
    :alt: BioConda Install
.. image:: https://img.shields.io/github/downloads/cpockrandt/phylocsfpp/total.svg
    :target: https://github.com/cpockrandt/phylocsfpp/releases/latest
    :alt: Github All Releases
.. image:: https://travis-ci.org/cpockrandt/phylocsfpp.svg?branch=master
    :target: https://travis-ci.org/cpockrandt/phylocsfpp
    :alt: Travis CI
.. image:: https://img.shields.io/badge/License-BSD%203--Clause-blue.svg
    :target: https://opensource.org/licenses/BSD-3-Clause
    :alt: BSD3 License

.. contents::
   :local:
   :depth: 2

Introduction
^^^^^^^^^^^^

PhyloCSF can identify protein-coding regions in the genome based on multiple-sequence alignments.
PhyloCSF++ is an implementation of the original methods `PhyloCSF`_ and `PhyloCSF HMM`_.
It allows you to easily create browser tracks for any genome to identify coding regions or score local alignments.

If you find our implementation useful, please consider citing it:

Christopher Pockrandt, Martin Steinegger, Steven Salzberg. **Title**. `Bioinformatics`_, 2021.

Please also consider citing the original developers of the method:

Lin MF, Jungreis I, and Kellis M. PhyloCSF: a comparative genomics method to distinguish protein-coding and non-coding regions. Bioinformatics, 2011.

.. _Bioinformatics: https://doi.org/10.1093/bioinformatics/btaa222

Installation
^^^^^^^^^^^^

Binaries
""""""""

.. Source of linux.svg: https://svgsilh.com/image/2025536.html

.. |VERSION| replace:: 1.3.0
.. |BUILD_DATE| replace:: 2020-06-17

Bioconda
""""""""

::

    $ conda install -c bioconda phylocsfpp

Building from source
""""""""""""""""""""

If you want to build it from source, we recommend cloning the git repository as shown below.
Please note that building from source can easily take 10 minutes and longer depending on your machine and compiler.

::

    $ git clone https://github.com/cpockrandt/PhyloCSFpp.git
    $ mkdir PhyloCSFpp/build && cd PhyloCSFpp/build
    $ cmake .. -DCMAKE_BUILD_TYPE=Release
    $ make

You can then install PhyloCSF++ as follows

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
  GNU scientific library (gsl)

Getting started
^^^^^^^^^^^^^^^

Building tracks
"""""""""""""""

TODO

::

    $ phylocsf++ ...

Scoring alignments
""""""""""""""""""

TODO

::

    $ phylocsf++ ...

Training my own model
"""""""""""""""""""""

TODO

::

    $ phylocsf++ ...

Help pages and examples
"""""""""""""""""""""""

A detailed list of arguments and explanations can be retrieved with ``--help``:

::

    $ phylocsf++ --help
    $ phylocsf++ tracks --help
    $ phylocsf++ regions --help

More detailed examples can be found in the `Wiki <https://github.com/cpockrandt/PhyloCSFpp/wiki>`_.

LICENSE
^^^^^^^

This is an implementation of the original methods (`PhyloCSF`_ and `PhyloCSF HMM`_), which were released under the GNU AGPL v3 and Apache License v2.
The algorithms were not changed except for running time improvements and were explicitly stated in the code.

TODO

.. _PhyloCSF: https://github.com/mlin/PhyloCSF
.. _PhyloCSF HMM: https://github.com/iljungr/PhyloCSFCandidateCodingRegions