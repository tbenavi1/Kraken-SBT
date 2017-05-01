from BitVector import BitVector
from collections import defaultdict
from ete3 import NCBITaxa
import mmh3
import operator
import os
import sys

class BloomFilter:
	
	def __init__(self, size, num_hashes, bitvector = None):
		self.size = size
		self.num_hashes = num_hashes
		if bitvector is None:
			self.bv = BitVector(intVal = 0, size = size)
		else:
			self.bv = bitvector
	
	def add(self, string):
		for seed in range(self.num_hashes):
			result = mmh3.hash128(string, seed) % self.size
			self.bv[result] = 1
	
	def contains(self, string):
		for seed in range(self.num_hashes):
			result = mmh3.hash128(string, seed) % self.size
			if self.bv[result] == 0:
				return False
		return True

#recall: command to generate name_ftpdirpaths from assembly_summary.txt (ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt) is given below:
#awk -F "\t" '$12=="Complete Genome" && $11=="latest"{print $8, $20}' assembly_summary.txt > name_ftpdirpaths
#each line of name_ftpdirpaths is space delimited; everything before the last space is the taxon name ($8), everything after the last space is the ftpdirpath ($20)
def get_taxonid_to_dumpsfilenames(name_ftpdirpaths_filename): 
	#search through all the data files and create a dictionary that maps taxonids to dumpsfilenames (dumps files are generated from Jellyfish and contain kmers)
	taxonid_to_dumpsfilenames = defaultdict(list) #a given taxonid may map to multiple dumpsfiles, thus each value in the dictionary is a list
	for line in open(name_ftpdirpaths_filename):
		splits = line.strip().split(' ') #split on each of the spaces
		name = ' '.join(splits[:-1]) #concatenate everything before the last space
		ftpdirpath = splits[-1] #everything after the last space
		dumpsfilename = './Bacteria_Genomes/' + ftpdirpath.split('/')[-1] + '_genomic.fna.gz.dumps'
		name_to_listtaxonid = ncbi.get_name_translator([name]) #a dictionary with name as key and [taxonid] as value
		listtaxonid = [taxonid for [taxonid] in name_to_listtaxonid.values()]
		if listtaxonid == []:
			if name == 'Donghicola sp. JLT3646': #upon further inspection, this name in assembly_summary.txt has been updated in NCBI
				print('Changing ' + name + ' to Marivivens sp. JLT3646')
				sys.stdout.flush()
				name = 'Marivivens sp. JLT3646'
				name_to_listtaxonid = ncbi.get_name_translator([name])
				listtaxonid = [taxonid for [taxonid] in name_to_listtaxonid.values()]
				taxonid = listtaxonid[0]
				taxonid_to_dumpsfilenames[taxonid].append(dumpsfilename)
			if name == 'Mycobacterium intracellulare MOTT-64': #upon further inspection, this name in assembly_summary.txt is curiously absent from NCBI
				print(name+' is not in the NCBI database.')
				sys.stdout.flush()
		else:
			taxonid = listtaxonid[0]
			taxonid_to_dumpsfilenames[taxonid].append(dumpsfilename)
	return taxonid_to_dumpsfilenames

def bf_from_bvfilename(bvfilename):
	bitvector = BitVector(filename = bvfilename) #loads file but doesn't actually read the data
	numbits = os.stat(bvfilename).st_size * 8
	bitvector = bitvector.read_bits_from_file(numbits) #reads the data into the variable bitvector
	bf = BloomFilter(len(bitvector), 3, bitvector) #We use num_hashes = 3
	return bf

def get_tree(name_ftpdirpaths_filename,num_taxons = 0):
	ncbi = NCBITaxa()
	taxonid_to_dumpsfilenames = get_taxonid_to_dumpsfilenames(name_ftpdirpaths_filename)
	
	#get the desired number of unique taxonids, in order to create the phylogeny tree
	taxonids = taxonid_to_dumpsfilenames.keys() #length is 6,740 as expected (6,741 minus the one not in NCBI)
	taxonids = list(set(taxonids)) #4,526 taxonids are unique
	if num_taxons != 0:
		taxonids = taxonids[:num_taxons] #smaller set of taxonids for tree construction and testing
	
	#return desired phylogeny tree
	return ncbi.get_topology(taxonids) #5,360 total nodes

if __name__=="__main__":
	
	tree = get_tree('name_ftpdirpaths')
	tree_test = get_tree('name_ftpdirpaths', 10)
	
	#construct the bloomfilters (only necessary for the first time building the database)
	#actually, the end user never needs to perform this step, since they will download the bloom filters from the beginning
	#construct_bloomfilters(tree)
	#construct_bloomfilters(tree_test)
