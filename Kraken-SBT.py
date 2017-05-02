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
	taxonid_to_dumpsfilenames = get_taxonid_to_dumpsfilenames(name_ftpdirpaths_filename)
	
	#get the desired number of unique taxonids, in order to create the phylogeny tree
	taxonids = taxonid_to_dumpsfilenames.keys() #length is 6,740 as expected (6,741 minus the one not in NCBI)
	taxonids = list(set(taxonids)) #4,526 taxonids are unique
	if num_taxons != 0:
		taxonids = taxonids[:num_taxons] #smaller set of taxonids for tree construction and testing
	
	#return desired phylogeny tree
	return taxonid_to_dumpsfilenames, ncbi.get_topology(taxonids) #5,360 total nodes for full dataset

def construct_bloomfilters(tree):
		
	edited_names = []
	i=0
	numnodes = len(list(tree.traverse()))
	for node in tree.traverse('postorder'): #traverse each node in the tree from leaf level upward to root level
		i+=1
		taxonid = int(node.name)
		name = ncbi.translate_to_names([taxonid])[0]
		print('Node ' + str(i) + ' out of ' + str(numnodes) + ': ' + name)
		sys.stdout.flush()
		if "/" in name:
			edited_names.append(name)
		edited_name = name.replace(' ', '_').replace('/', '_') #replace spaces and / with underscores so that it can be a filename, e.g. Acinetobacter calcoaceticus/baumannii complex
		bv_filename = edited_name + '.bv'
		if False: #os.path.exists(bv_filename):
			print(bv_filename + ' already exists. Loading...')
			sys.stdout.flush()
			node.bf = bf_from_bvfilename(bv_filename)
		else:
			node.kmers = []
			if taxonid in taxonid_to_dumpsfilenames:
				dumpsfilenames = taxonid_to_dumpsfilenames[taxonid]
				numfiles = len(dumpsfilenames)
				j = 0
				for dumpsfilename in dumpsfilenames:
					j+= 1
					print('Node ' + str(i) + ' out of ' + str(numnodes) + ': ' + name + '; File ' + str(j) + ' out of ' + str(numfiles))
					sys.stdout.flush()
					for line in open(dumpsfilename):
						kmer = line.strip().split(' ')[0]
						node.kmers.append(kmer)
			if not node.is_leaf():
				children = node.children
				for child in children:
					node.kmers.extend(child.kmers)
					delattr(child, 'kmers')
			size = int(round(len(node.kmers) * 4.25,-3))
			node.bf = BloomFilter(size, num_hashes = 3)
			for kmer in node.kmers:
				node.bf.add(kmer)
			f = open(bv_filename, 'wb')
			node.bf.bv.write_to_file(f)
			f.close()
			delattr(node, 'bf')

#query the tree
def query_tree(tree, querytaxonid):
	
	def get_query_kmers(querytaxonid):
		queryname = ncbi.translate_to_names([querytaxonid])[0]
		print('Query name is ' + queryname)
		sys.stdout.flush()
		querydumpsfilenames = taxonid_to_dumpsfilenames[querytaxonid]
		querykmers = []
		for querydumpsfilename in querydumpsfilenames:
			for line in open(querydumpsfilename):
				kmer = line.strip().split(' ')[0]
				querykmers.append(kmer)
		return queryname, querykmers
	
	def get_next_node_kmers(children, current_kmers, threshold):
		node_kmers = []
		for child in children:
			taxonid = int(child.name)
			name = ncbi.translate_to_names([taxonid])[0]
			edited_name = name.replace(' ', '_').replace('/', '_')
			bv_filename = edited_name + '.bv'
			print('Loading ' + name)
			sys.stdout.flush()
			child.bf = bf_from_bvfilename(bv_filename)
			print(name + ' loaded')
			sys.stdout.flush()
			kmer_matches = []
			for kmer in current_kmers:
				if child.bf.contains(kmer):
					kmer_matches.append(kmer)
			delattr(child, 'bf')
			if len(kmer_matches) > threshold:
				node_kmers.append((child, kmer_matches))
		return node_kmers
	
	responses = {}
	queryname, querykmers = get_query_kmers(querytaxonid)
	num_kmers = len(querykmers)
	threshold = num_kmers * 0.5
	node_kmers_to_query = []
	node_kmers_to_query.append((tree,querykmers))
	while node_kmers_to_query:
		current_node, current_kmers = node_kmers_to_query.pop()
		current_taxonid = int(current_node.name)
		current_name = ncbi.translate_to_names([current_taxonid])[0]
		print('Current node: ' + current_name)
		sys.stdout.flush()
		if current_node.is_leaf():
			responses[current_name] = len(current_kmers)/num_kmers
			print('Proportion of query kmers matching ' + current_name + ': ' + str(len(current_kmers)/num_kmers))
			sys.stdout.flush()
		else:
			children = current_node.children
			node_kmers = get_next_node_kmers(children, current_kmers, threshold)
			if node_kmers == []:
				responses[current_name] = len(current_kmers)/num_kmers
				print('Proportion of query kmers matching ' + current_name + ': ' + str(len(current_kmers)/num_kmers))
				sys.stdout.flush()
			else:
				node_kmers_to_query.extend(node_kmers)
	return sorted(responses.items(), key = operator.itemgetter(1), reverse = True)

if __name__=="__main__":
	
	ncbi = NCBITaxa()
	
	num_taxons = sys.argv[1]
	
	taxonid_to_dumpsfilenames, tree = get_tree('name_ftpdirpaths', num_taxons)
	
	command = sys.argv[2]
	
	if command == "bloomfilter":
		#construct the bloomfilters (only necessary for the first time building the database)
		#actually, the end user never needs to perform this step, since they will download the bloom filters from the beginning
		construct_bloomfilters(tree)
	
	if command == "query":
		querytaxonid = sys.argv[3]
		print(query_tree(tree, querytaxonid))
