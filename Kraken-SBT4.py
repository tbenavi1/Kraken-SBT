from BitVector import BitVector
from collections import defaultdict
from ete3 import NCBITaxa
import mmh3
import numpy as np
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
		dumpsfilename = './Bacteria_Genomes/' + ftpdirpath.split('/')[-1] + '_genomic.fna.dumps'
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

def add_bloomfilters(tree):
		
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
def query_tree(querytaxonid, tree):
	
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
	taxonid_to_dumpsfilenames = get_taxonid_to_dumpsfilenames('name_ftpdirpaths')
	
	#get the list of unique taxonids, in order to create the phylogeny tree
	taxonids = taxonid_to_dumpsfilenames.keys() #length is 6,740 as expected (6,741 minus the one not in NCBI)
	unique_taxonids = list(set(taxonids)) #4,526 taxonids are unique
	taxonids_test = unique_taxonids[:10] #smaller set of taxonids for tree construction and testing

	#create full phylogeny tree and test phylogeny tree
	tree = ncbi.get_topology(unique_taxonids) #5,360 total nodes
	tree_test = ncbi.get_topology(taxonids_test)

	#construct (or load if already exist) the bloomfilters
	add_bloomfilters(tree)
	#add_bloomfilters(tree_test)

def get_Kraken_accuracy(fastafilename,kraken_response_taxonid):
	kraken_response_name = ncbi.translate_to_names([kraken_response_taxonid])
	true_response_taxonid = [taxonid for taxonid in taxonid_to_dumpsfilenames if fastafilename + '.dumps' in taxonid_to_dumpsfilenames[taxonid]]
	true_response_name = ncbi.translate_to_names(true_response_taxonid)
	return 'Fasta filename: ' + fastafilename + '; Kraken response taxonid: ' + str(kraken_response_taxonid) + '; Kraken response name: ' + str(kraken_response_name) + '; true response taxonid; ' + str(true_response_taxonid) + '; true response name: ' + str(true_response_name)


#Kraken needs 160GB of disk space, and 75GB of RAM (for bacterial, viral, and )
#Bacteria_Genomes/GCF_900169565.1_NSJP_Ch1_genomic.fna 330214 = Nitrospira defulvii; actually 1325564 = Nitrospira japonica
#"Fasta filename: ./Bacteria_Genomes/GCF_900170005.1_PRJEB19975_genomic.fna; Kraken response taxonid: 657308; Kraken response name: ['Gordonibacter pamelaeae 7-10-1-b']; true response taxonid; [1841863]; true response name: ['Gordonibacter sp. Marseille-P2775']"

#2 sequences
#"Fasta filename: ./Bacteria_Genomes/GCF_900169085.1_SWA-2_genomic.fna; Kraken response taxonid: 813; Kraken response name: ['Chlamydia trachomatis']; true response taxonid; [83559]; true response name: ['Chlamydia suis']"
#"Fasta filename: ./Bacteria_Genomes/GCF_900169085.1_SWA-2_genomic.fna; Kraken response taxonid: 707187; Kraken response name: ['Chlamydia trachomatis G/11074']; true response taxonid; [83559]; true response name: ['Chlamydia suis']"

#"Fasta filename: ./Bacteria_Genomes/GCF_000266945.1_ASM26694v1_genomic.fna; Kraken response taxonid: 706587; Kraken response name: ['Desulfomonile tiedjei DSM 6799']; true response taxonid; [706587]; true response name: ['Desulfomonile tiedjei DSM 6799']"
#"Fasta filename: ./Bacteria_Genomes/GCF_000952915.1_9231-Abomsa_assembly_V8_genomic.fna; Kraken response taxonid: 340047; Kraken response name: ['Mycoplasma capricolum subsp. capricolum ATCC 27343']; true response taxonid; [40480]; true response name: ['Mycoplasma capricolum subsp. capripneumoniae']"
#Bacteria_Genomes/GCF_001511815.1_TUE45_genomic.fna Streptomyces hygroscopicus subsp. jinggangensis 5008 true response taxonid; [1926]; true response name: ['Streptomyces reticuli']"
#GCF_001457595.1_NCTC8159_genomic.fna Mycobacterium smegmatis str. MC2 155 true response taxonid; [1772]; true response name: ['Mycobacterium smegmatis']"




#3 had wrong identification
#[ncbi.translate_to_names([x]) for x in ncbi.get_lineage([x for [x] in ncbi.get_name_translator(['Buchnera aphidicola']).values()][0])]

#Acinetobacter_baumannii_ZW85-1.bf needs to be redone
#edited_names = ['Acinetobacter calcoaceticus/baumannii complex', 'Neisseria meningitidis H44/76', 'Rhizobium/Agrobacterium group', 'Listeria monocytogenes serotype 1/2a str. 88-0478', 'Listeria monocytogenes serotype 1/2a str. 08-6569', 'Listeria monocytogenes serotype 1/2a str. 08-6997', 'Listeria monocytogenes serotype 1/2a str. 10-0815', 'Listeria monocytogenes serotype 1/2a str. 10-1047', 'Listeria monocytogenes serotype 1/2a']


#Most populous in database
#Bordella pertussis 286
#Escherichia coli 188
#Staphylococcus aureus 73
#Campylobacter jejuni 72
#Klebsiella pneuniae 70
#Pseudomonas aeruginosa 60
#Neisseria meningitidis 59
#Listeria monocytogenes 57
#Corynebacterium pseudotuberculosis 44

		#matches = [sum([child.bf.__contains__(kmer) for kmer in querykmers]) for child in children]
		#print(matches)
		#new_current_node_index = np.argmax(matches)
		#current_node = current_node.children[new_current_node_index]
#1 hr for kraken to do 800000 sequences from i100 dataset
#2 hr to classify everything
