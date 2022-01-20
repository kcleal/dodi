====
dodi
====

dodi chooses an optimal set of alignments from a candidate list. This can be used for
improving split-read mapping, or for selectively mapping to target regions of the genome

Installation
------------
Install using::

    $ git clone https://github.com/kcleal/dodi.git; cd dodi
    $ python setup.py install

To satisfy requirements::

    $ pip install -r requirements.txt; pip install .

Run tests::

    $ dodi test

Requires Python>=3.6, cython and >=c++11 compiler.
Python packages needed are listed in requirements.txt.


Usage
-----
Given a list of candidate alignments in `.sam` format, `dodi` will select an optimal spanning set of
alignments. Common usage is to run `dodi` in a pipe, downstream of a mapper such as bwa mem or Last.

To use with bwa mem, make sure the `-a` flag is used, to report all alignments::

    $ bwa mem -a -t8 ref.fa read1.fq read2.fq | dodi - > output_alignments.sam


