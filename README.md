# Kraken-SBT

## Prerequisites
1. Python 3 from https://www.python.org/downloads/release/python-361/
2. ete3 (for NCBI taxonomy tree) from http://etetoolkit.org/download/
3. bitarray from https://pypi.python.org/pypi/bitarray/
4. mmh3 from https://pypi.python.org/pypi/mmh3/2.3.1
5. Jellyfish 2 from http://www.genome.umd.edu/jellyfish.html

## Setting up your environment
Running:
```bash
sh setup.sh
```
will first download the 6,741 complete RefSeq bacterial genomes into a directory called Bacteria_Genomes. Then, it will run jellyfish count and jellyfish dump to get the kmer dumps files for each of the genomes.
Note: This will take a very long time and it is not necessary for the typical end user. This is only necessary for the developer.

## Running for the first time
Kraken-SBT can run on the full dataset (~4500 taxonids), or on a subset of the taxonids. To build the tree on the full dataset, run:
```bash
python Kraken-SBT buildtree 0
```
To build the tree on a subset, the first argument of buildtree should be the desired number of taxonids. For instance:
```bash
python Kraken-SBT buildtree 10
```
## Construct Bloom Filter Database Files
To construct the bloomfilters for the first time, run:
```bash
python Kraken-SBT bloomfilters 0
```
Running on the entire dataset will take an inordinate amount of time and memory. This is only necessary for the developer.

## Querying the tree
To query the full tree for a given taxonid (e.g. 385025), and return the nodes for which the proportion of matching query kmers exceeds a given threshold, run:
```bash
python Kraken-SBT query 0 385025 threshold
```
To query the full tree for a given jellyfish dumps file containing a set of kmers, and return the nodes for which the proportion of matching query kmers exceeds a given threshold, run:
```bash
python Kraken_SBT query 0 queryfilename threshold
```
