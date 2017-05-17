# Kraken-SBT

## Prerequisites
1. Python 3 from https://www.python.org/downloads/release/python-361/
2. ete3 (for NCBI taxonomy tree) from http://etetoolkit.org/download/
3. bitarray from https://pypi.python.org/pypi/bitarray/
4. mmh3 from https://pypi.python.org/pypi/mmh3/2.3.1
5. jellyfish 2 from http://www.genome.umd.edu/jellyfish.html

## Setting up your environment
Running:
```bash
sh setup.sh
```
will first download the 6,741 complete RefSeq bacterial genomes into a directory called Bacteria_Genomes. Then, it will run jellyfish count and jellyfish dump to get the kmer dumps files for each of the genomes.
Note: This will take a very long time and it is not necessary for the typical end user. This is only necessary for the developer.

## Running for the first time
Kraken-SBT can run on the tree built with the complete RefSeq genomes (~4500 taxonids), or on the tree built with the i100 metagenomics genomes (56 taxonids). Additionally, Kraken-SBT can run on the entire tree (all taxonids), or on a subset of the tree (a given number of taxonids). To build the tree with the complete RefSeq genomes on all taxonids, run:
```bash
python Kraken-SBT buildtree complete 0
```
To build the tree with the complete RefSeq genomes on a subset of the taxonids, the last argument should be the desired number of taxonids. For instance:
```bash
python Kraken-SBT buildtree complete 10
```
To build the tree with the i100 genomes on all taxonids, run:
```bash
python Kraken-SBT buildtree i100 0
```
To build the tree with the i100 genomes on a subset of the taxonids, the last argument should be the desired number of taxonids. For instance:
```bash
python Kraken-SBT buildtree i100 10
```
## Construct Bloom Filter Database Files
To construct the bloomfilters (only necessary to do once), you can use the same arguments as for buildtree above. For example, to construct the bloomfilters on first 10 taxonids of the complete RefSeq tree, run:
```bash
python Kraken-SBT bloomfilters complete 0
```
Running on the entire dataset will take an inordinate amount of time and memory. This is only necessary for the developer. Note: Subsequent calls to bloomfilters will recognize if the bloom filter already exists, and thus will not waste time recreating the bloom filter. To test the construction of the bloom filters, you can run on increasingly larger trees.

## Querying the tree
To query the complete tree on all taxonids for a given number of kmers (num_kmers) of a given taxonid (e.g. 385025), and return the nodes for which the proportion of matching query kmers exceeds a given threshold (theta), run:
```bash
python Kraken-SBT query complete 0 385025 num_kmers theta num_workers
```
The last argument (num_workers) is the number of processing threads used for the computation.
To query the complete tree on all taxonids for a given number of kmers of a given jellyfish dumps file containing a set of kmers, and return the nodes for which the proportion of matching query kmers exceeds a given threshold (theta), run:
```bash
python Kraken_SBT query complete 0 queryfilename num_kmers theta num_workers
```
The last argument is the number of processing threads used for the computation, as before.
