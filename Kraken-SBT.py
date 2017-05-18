from bitarray import bitarray
from collections import defaultdict
from ete3 import NCBITaxa
import mmh3
import multiprocessing as mp
import operator
import os
import queue
import sys
import threading

#define the BloomFilter class
class BloomFilter:
	
	def __init__(self, size, num_hashes, bitvector = None):
		self.size = size
		self.num_hashes = num_hashes
		if bitvector is None:
			self.bv = bitarray(size) #initialize a bitvector with the given size
			self.bv.setall(False) #initialize the bitvector to be all 0s
		else:
			self.bv = bitvector #if a bitvector is provided, use it
	
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
def get_taxonid_to_readfilenames(name_ftpdirpaths_filename): #searches through all the data files and create a dictionary that maps taxonids to readfilenames
	ncbi = NCBITaxa()
	taxonid_to_readfilenames = defaultdict(list) #a given taxonid may map to multiple readfiles, thus each value in the dictionary is a list
	for line in open(name_ftpdirpaths_filename):
		splits = line.strip().split(' ') #split on each of the spaces
		name = ' '.join(splits[:-1]) #concatenate everything before the last space
		ftpdirpath = splits[-1] #everything after the last space
		readfilename = './Bacteria_Genomes/' + ftpdirpath.split('/')[-1] + '_genomic.fna.gz'
		name_to_taxonid = ncbi.get_name_translator([name]) #a dictionary with name as key and [taxonid] as value
		listtaxonid = [taxonid for [taxonid] in name_to_taxonid.values()] #if the name was found in the ncbi database, listtaxonid is a list containing the taxonid; otherwise, an empty list.
		if listtaxonid:
			taxonid = listtaxonid[0]
			taxonid_to_readfilenames[taxonid].append(readfilename)
		else:
			if name == 'Donghicola sp. JLT3646': #upon further inspection, this name in assembly_summary.txt has been updated in NCBI
				print('Changing', name, 'to Marivivens sp. JLT3646')
				name = 'Marivivens sp. JLT3646'
				name_to_taxonid = ncbi.get_name_translator([name])
				listtaxonid = [taxonid for [taxonid] in name_to_taxonid.values()]
				taxonid = listtaxonid[0]
				taxonid_to_readfilenames[taxonid].append(readfilename)
			if name == 'Mycobacterium intracellulare MOTT-64': #upon further inspection, this name in assembly_summary.txt is curiously absent from NCBI
				print(name, 'is not in the NCBI database')
	return taxonid_to_readfilenames

def get_tree(taxonids, num_taxonids = 0):
	ncbi = NCBITaxa()
	
	if num_taxonids != 0:
		taxonids = taxonids[:num_taxonids] #smaller set of taxonids for tree construction and testing
	
	return ncbi.get_topology(taxonids) #5,360 total nodes for full dataset

def get_taxonid_to_name(tree):
	ncbi = NCBITaxa()
	taxonid_to_name = {}
	
	for node in tree.traverse():
		taxonid = int(node.name)
		taxonid_to_name[taxonid] = ncbi.translate_to_names([taxonid])[0]
	
	return taxonid_to_name

def write_descendantfiles(tree, taxonid_to_name):
	num_nodes = len(list(tree.traverse()))
	for i, node in enumerate(tree.traverse()):
		taxonid = int(node.name)
		name = taxonid_to_name[taxonid]
		print('Node', i+1, 'out of', str(num_nodes) + ':', name)
		edited_name = name.replace(' ', '_').replace('/', '_')  #replace spaces and / with underscores so that it can be a filename, e.g. Acinetobacter calcoaceticus/baumannii complex
		descendantfilenames_filename = edited_name + '.descendantfilenames'
		readfilenames = taxonid_to_readfilenames[taxonid]
		numfiles = len(readfilenames)
		for j, readfilename in enumerate(readfilenames):
			print('Node', i+1, 'out of', str(num_nodes) + ':', name + '; File', j+1, 'out of', numfiles)
			f = open(descendantfilenames_filename,'a')
			f.write(readfilename + '\n')
			f.close()
		if not node.is_leaf():
			num_descendants = len(node.get_descendants())
			for k, descendant in enumerate(node.get_descendants()):
				descendant_taxonid = int(descendant.name)
				descendant_name = taxonid_to_name[descendant_taxonid]
				print('Node', i+1, 'out of', str(num_nodes) + ':', name + '; Descendant Node', k+1, 'out of', str(num_descendants) + ':', descendant_name)
				descendant_filenames = taxonid_to_readfilenames[descendant_taxonid]
				descendant_numfiles = len(descendant_filenames)
				for l, descendant_filename in enumerate(descendant_filenames):
					print('Node', i+1, 'out of', str(num_nodes) + ':', name + '; Descendant Node', k+1, 'out of', str(num_descendants) + ':', descendant_name + '; Descendant File', l+1, 'out of', descendant_numfiles)
					f = open(descendantfilenames_filename,'a')
					f.write(descendant_filename + '\n')
					f.close()

def read_bloomfiltersizes(bloomfiltersizes_filename): #The bloomfiltersizes file was created by running the KMC kmer counter on the files of descendant filenames created from the function write_descendantfiles (above)
	bloomfiltersizes = {}
	for line in open(bloomfiltersizes_filename):
		edited_name, size = line.strip().split(' ')
		bloomfiltersizes[edited_name] = int(size)
	return bloomfiltersizes

def construct_bloomfilter(taxonid_to_name, bloomfiltersizes, num_nodes, nodesQueue):
	while True:
		i, node = nodesQueue.get()
		if node is None:
			break
		taxonid = int(node.name)
		name = taxonid_to_name[taxonid]
		print('Constructing bloom filter for node', i+1, 'out of', str(num_nodes) + ':', name)
		edited_name = name.replace(' ', '_').replace('/', '_')
		descendantfilenames_filename = edited_name + '.descendantfilenames'
		bv_filename = './bloomfilters/' + edited_name + '.bv'
		if os.path.exists(bv_filename):
			print(bv_filename, 'already exists.')
			nodesQueue.task_done()
		else:
			size = bloomfiltersizes[edited_name] #the size was determined by counting unique kmers using kmc kmer counter
			size = int(round(size * 4.25,-3)) #multiply by 4.25; round to nearest 1000 (for a 13% false positive rate; bitvector length must be a multiple of 8)
			node.bf = BloomFilter(size, num_hashes = 3) #num_hashes = 3 (for a 13% positive rate)
			with open(descendantfilenames_filename) as f:
				num_files = sum(1 for _ in f)
			for j, descendant_filename in enumerate(open(descendantfilenames_filename)): #this for loop and the for loop below add all the kmers from all the descendant node files to the bloom filter
				descendant_filename = descendant_filename.strip() + '.dumps'
				print('Reading file', j+1, 'out of', str(num_files) + ':', name)
				for line in open(descendant_filename):
					kmer = line.strip().split(' ')[0]
					node.bf.add(kmer)
			f = open(bv_filename, 'wb') #open bitvector file for writing
			node.bf.bv.tofile(f) #write bitvector to file
			f.close()
			print('Bloom filter for', name, 'constructed.')
			delattr(node, 'bf') #remove bloomfilter from memory
			nodesQueue.task_done()

def construct_bloomfilters(tree, taxonid_to_name, bloomfiltersizes):
	num_nodes = len(list(tree.get_descendants()))
	
	nodesManager = mp.Manager()
	nodesQueue = nodesManager.Queue()
	
	num_processes = 100
	processes = []
	for i in range(num_processes):
		process = mp.Process(target = construct_bloomfilter, args = (taxonid_to_name, bloomfiltersizes, num_nodes, nodesQueue))
		process.start()
		processes.append(process)
	
	for i, node in enumerate(tree.get_descendants()): #don't need to get a bloom filter for the root node 'Bacteria'
		nodesQueue.put((i,node))
	
	nodesQueue.join()
	
	for i in range(num_processes):
		nodesQueue.put((None, None))
	for p in processes:
		p.join()

def get_query_kmers(taxonid_to_readfilenames, taxonid_to_name, query, max_num_kmers, num_workers):
	querykmers = []
	try: #if the query is a taxonid
		querytaxonid = int(query)
		queryname = taxonid_to_name[querytaxonid]
		print('Query name is', queryname)
		queryreadfilenames = taxonid_to_readfilenames[querytaxonid]
		for queryreadfilename in queryreadfilenames:
			querydumpsfilename = queryreadfilename + '.dumps'
			with open(querydumpsfilename) as f:
				for line in f:
					if i == max_num_kmers and i != 0: #only get at most max_num_kmers kmers of the query; if max_num_kmers == 0, get all kmers of the query
						break
					kmer = line.strip().split(' ')[0]
					querykmers.append(kmer)
	except: #if the query is a filename
		queryfilename = query
		print('Query filename is', queryfilename)
		with open(queryfilename) as f:
			for i, line in enumerate(f):
				if i == max_num_kmers and i != 0: #only get at most max_num_kmers kmers of the query; if max_num_kmers == 0, get all kmers of the query
					break
				kmer = line.strip().split(' ')[0]
				querykmers.append(kmer)
	return querykmers

def bf_from_bvfilename(bv_filename):
	f = open(bv_filename, 'rb') #open file with bitvector
	bitvector = bitarray() #initalize bitvector
	bitvector.fromfile(f) #load bitvector from file
	f.close()
	bf = BloomFilter(bitvector.length(), 3, bitvector) #create bloom filter with num_hashes = 3
	return bf

def get_next_node_kmers(taxonid_to_name, threshold_proportion, num_workers, num_kmers, current_kmers, p, children):
	next_node_kmers = [] #a list of (node, kmers) tuples
	
	#Load the bloom filters
	for i, child in enumerate(children):
		taxonid = int(child.name)
		name = taxonid_to_name[taxonid]
		edited_name = name.replace(' ', '_').replace('/', '_')
		bv_filename = './bloomfilters/' + edited_name + '.bv'
		print('Loading', name)
		children[i].bf = bf_from_bvfilename(bv_filename)
		print(name, 'loaded')
	
	kmer_matches = [[] for _ in range(len(children))]
	f = 0.13 #the false positive rate of the bloom filters	

	#add matching kmers if exceed threshold_proportion
	for i, child in enumerate(children):
		for current_kmer in current_kmers:
			if child.bf.contains(current_kmer):
				kmer_matches[i].append(current_kmer)
		q = len(kmer_matches[i])/num_kmers
		adjusted_proportion = (q-p*f)/(1-f)
		if adjusted_proportion >= threshold_proportion:
			next_node_kmers.append((child, kmer_matches[i], adjusted_proportion))
	
	#remove bloom filters from memory
	for i, child in enumerate(children):
		delattr(child, 'bf')
	
	return next_node_kmers

def analyze_node(taxonid_to_name, threshold_proportion, num_workers, name_to_proportion, num_kmers, workQueue, queriedQueue): #each process targets this function to analyze its current node
	while True:
		current_node, current_kmers, adjusted_proportion = workQueue.get()
		if current_node is None:
			break
		p = len(current_kmers)/num_kmers
		current_taxonid = int(current_node.name)
		current_name = taxonid_to_name[current_taxonid]
		print('Current node:', current_name)
		if current_node.is_leaf():
			name_to_proportion[current_name] = adjusted_proportion
			print('Proportion of query kmers matching', current_name + ':', adjusted_proportion)
			workQueue.task_done()
		else:
			children = current_node.children
			next_node_kmers = get_next_node_kmers(taxonid_to_name, threshold_proportion, num_workers, num_kmers, current_kmers, p, children)
			num_queried = queriedQueue.get()
			num_queried += len(children)
			queriedQueue.put(num_queried)
			if next_node_kmers == []: # if all of the children are a dead end (i.e. proportion matching doesn't exceed threshold)
				name_to_proportion[current_name] = adjusted_proportion
				print('Proportion of query kmers matching', current_name + ':', adjusted_proportion)
				workQueue.task_done()
			else: #if at least one child is not a dead end
				for (child, kmer_matches, adjusted_proportion) in next_node_kmers: #add matching children to work queue
					workQueue.put((child, kmer_matches, adjusted_proportion))
					print('Added to queue.')
				workQueue.task_done()

#query the tree
def query_tree(taxonid_to_readfilenames, tree, taxonid_to_name, query, max_num_kmers, threshold_proportion, num_workers):
	
	workManager = mp.Manager() #a multiprocessing manager to manage all the shared data between workers
	
	name_to_proportion = workManager.dict() # a multiprocessing dictionary with name as key and the proportion of kmers matching the name as the value; since it is multiprocessing, it can be accessed and modified by all workers
	
	#get query kmers
	querykmers = get_query_kmers(taxonid_to_readfilenames, taxonid_to_name, query, max_num_kmers, num_workers)
	num_kmers = len(querykmers)
	print('Number of kmers in query:', num_kmers)
	
	#get the total number of nodes in the tree
	num_nodes = len(list(tree.traverse()))
	
	#create and initialize the work queue
	workQueue = workManager.Queue() #this multiprocessing 'work' queue contains all the nodes and corresponding matching kmers for each query path down the tree
	workQueue.put((tree,querykmers, 1)) #'tree' corresponds to the root node of the tree; 1 refers to the (adjusted) proportion of kmers that match at this node
	
	queriedQueue = workManager.Queue() #this multiprocessing queue contains an integer corresponding to the number of nodes that have been visited during the querying process (summed along all paths)
	queriedQueue.put(0)
	
	#create new processes
	processes = []
	for _ in range(num_workers): #start processes to target 'analyze_node' function
		process = mp.Process(target = analyze_node, args = (taxonid_to_name, threshold_proportion, num_workers, name_to_proportion, num_kmers, workQueue, queriedQueue))
		process.start()
		processes.append(process)
	
	#wait for the work queue to be complete
	workQueue.join()
	
	#Add sentinel (stop) values to work queue to signal processes to stop, then wait for all processes to complete
	for i in range(num_workers):
		workQueue.put((None, None, None))
	for p in processes:
		p.join()
	
	num_queried = queriedQueue.get()
	print(num_queried, 'nodes queried out of', num_nodes)
	return sorted(name_to_proportion.items(), key = operator.itemgetter(1), reverse = True)

if __name__=="__main__":
	
	try:
		command = sys.argv[1]
	except:
		command = []
		print('Possible commands are \'buildtree\', \'bloomfilters\', and \'query\'.')
	
	if command:
		tree_to_build = sys.argv[2]
		num_taxonids = int(sys.argv[3])
		taxonid_to_readfilenames = get_taxonid_to_readfilenames('name_ftpdirpaths')
		
		#get the unique taxonids, in order to create the phylogeny tree
		taxonids = taxonid_to_readfilenames.keys() #length is 6,740 as expected (6,741 minus the one not in NCBI)
		taxonids = list(set(taxonids)) #4,526 taxonids are unique
		
		if tree_to_build == 'complete':
			tree = get_tree(taxonids, num_taxonids)
		
		if tree_to_build == 'i100':
			#tree built on i100 species
			#Note: Of the 85 species present in the i100 metagenomics sample, 3 are not bacteria, and 26 are not labeled as "Complete Genome" in RefSeq. Thus, this tree consists of the 56 complete bacteria species.
			tree = get_tree([470, 1392, 1396, 79880, 1423, 1428, 818, 518, 9, 28450, 199, 83558, 813, 1520, 1491, 43771, 777, 106590, 43989, 562, 263, 35554, 727, 731, 135577, 1584, 1624, 1358, 29546, 174, 1764, 1765, 1781, 164757, 1773, 485, 487, 529, 1219, 303, 316, 1076, 782, 1515, 28901, 60481, 382, 1280, 1311, 1313, 1314, 32046, 46541, 274, 51229, 632], num_taxonids)
		
		num_nodes = len(list(tree.traverse()))
		print('Tree has', num_nodes, 'nodes.')
		
		taxonid_to_name = get_taxonid_to_name(tree)
	
	if command == "descendantfilenames":
		f = open('editednames'):
		for node in tree.traverse():
			taxonid = int(node.name)
			name = taxonid_to_name[taxonid]
			edited_name = name.replace(' ', '_').replace('/', '_')
			f.write(edited_name + '\n')
		f.close()
		write_descendantfiles(tree, taxonid_to_name)
	
	if command == "bloomfilters":
		#construct the bloomfilters (only necessary for the first time building the database)
		#actually, the end user never needs to perform this step, since they will download the bloom filters directly
		bloomfiltersizes = read_bloomfiltersizes('bloomfiltersizes')
		construct_bloomfilters(tree, taxonid_to_name, bloomfiltersizes)
	
	if command == "query":
		query = sys.argv[4]
		max_num_kmers = int(sys.argv[5])
		threshold_proportion = float(sys.argv[6])
		num_workers = int(sys.argv[7])
		matches = query_tree(taxonid_to_readfilenames, tree, taxonid_to_name, query, max_num_kmers, threshold_proportion, num_workers)
		print(matches)
