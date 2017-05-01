# Kraken-SBT

## Prerequisites
1. Python 3 from https://www.python.org/downloads/release/python-361/
2. ete3 (for NCBI taxonomy tree) from http://etetoolkit.org/download/
3. BitVector from https://pypi.python.org/pypi/BitVector/3.4.7
4. mmh3 from https://pypi.python.org/pypi/mmh3/2.3.1
5. Jellyfish 2 from http://www.genome.umd.edu/jellyfish.html

## Setting up your environment
Running:
~~~bash
sh setup.sh
~~~
will first download the 6,741 complete RefSeq bacterial genomes into a directory called Bacteria_Genomes. Then, it will run jellyfish count and jellyfish dump to get the kmer dumps files for each of the genomes.
Note: This will take a very long time and it is not necessary for the typical end user. This is only necessary for the developer.

## Running for the first time
