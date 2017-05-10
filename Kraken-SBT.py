from bitarray import bitarray
from collections import defaultdict
from ete3 import NCBITaxa
from multiprocessing import Pool
import mmh3
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
		#readfilename = './Bacteria_Genomes/' + ftpdirpath.split('/')[-1] + '_genomic.fna.gz'
		readfilename = './Bacteria_Genomes/' + ftpdirpath.split('/')[-1] + '_genomic.fna'
		name_to_taxonid = ncbi.get_name_translator([name]) #a dictionary with name as key and [taxonid] as value
		listtaxonid = [taxonid for [taxonid] in name_to_taxonid.values()] #if the name was found in the ncbi database, listtaxonid is a list containing the taxonid; otherwise, an empty list.
		if listtaxonid:
			taxonid = listtaxonid[0]
			taxonid_to_readfilenames[taxonid].append(readfilename)
		else:
			if name == 'Donghicola sp. JLT3646': #upon further inspection, this name in assembly_summary.txt has been updated in NCBI
				print('Changing ' + name + ' to Marivivens sp. JLT3646')
				name = 'Marivivens sp. JLT3646'
				name_to_taxonid = ncbi.get_name_translator([name])
				listtaxonid = [taxonid for [taxonid] in name_to_taxonid.values()]
				taxonid = listtaxonid[0]
				taxonid_to_readfilenames[taxonid].append(readfilename)
			if name == 'Mycobacterium intracellulare MOTT-64': #upon further inspection, this name in assembly_summary.txt is curiously absent from NCBI
				print(name+' is not in the NCBI database')
	return taxonid_to_readfilenames

def bf_from_bvfilename(bv_filename):
	f = open(bv_filename, 'rb') #open file with bitvector
	bitvector = bitarray() #initalize bitvector
	bitvector.fromfile(f) #load bitvector from file
	f.close()
	bf = BloomFilter(bitvector.length(), 3, bitvector) #create bloom filter with num_hashes = 3
	return bf

def get_tree(taxonid_to_readfilenames,num_taxonids = 0):
	ncbi = NCBITaxa()
	
	#get the desired number of unique taxonids, in order to create the phylogeny tree
	taxonids = taxonid_to_readfilenames.keys() #length is 6,740 as expected (6,741 minus the one not in NCBI)
	taxonids = list(set(taxonids)) #4,526 taxonids are unique
	if num_taxonids != 0:
		taxonids = taxonids[:num_taxonids] #smaller set of taxonids for tree construction and testing
	
	#return desired phylogeny tree
	return ncbi.get_topology(taxonids) #5,360 total nodes for full dataset

def write_descendantfiles(tree):
	ncbi = NCBITaxa()
	num_nodes = len(list(tree.traverse()))
	for i, node in enumerate(tree.traverse()):
		taxonid = int(node.name)
		name = ncbi.translate_to_names([taxonid])[0]
		print('Node ' + str(i) + ' out of ' + str(num_nodes) + ': ' + name)
		edited_name = name.replace(' ', '_').replace('/', '_')  #replace spaces and / with underscores so that it can be a filename, e.g. Acinetobacter calcoaceticus/baumannii complex
		f = open("editednames", 'a')
		f.write(edited_name + '\n')
		f.close()
		descendantfilenames_filename = edited_name + '.descendantfilenames'
		if taxonid in taxonid_to_readfilenames:
			readfilenames = taxonid_to_readfilenames[taxonid]
			numfiles = len(readfilenames)
			j = 0
			for readfilename in readfilenames:
				j+= 1
				print('Node ' + str(i) + ' out of ' + str(num_nodes) + ': ' + name + '; File ' + str(j) + ' out of ' + str(numfiles))
				f = open(descendantfilenames_filename,'a')
				f.write(readfilename + '\n')
				f.close()
		if not node.is_leaf():
			num_descendants = len(node.get_descendants())
			k = 0
			for descendant in node.get_descendants():
				k += 1
				descendant_taxonid = int(descendant.name)
				descendant_name = ncbi.translate_to_names([descendant_taxonid])[0]
				print('Node ' + str(i) + ' out of ' + str(num_nodes) + ': ' + name + '; Descendant Node ' + str(k) + ' out of ' + str(num_descendants) + ': ' + descendant_name)
				descendant_filenames = taxonid_to_readfilenames[descendant_taxonid]
				descendant_numfiles = len(descendant_filenames)
				l = 0
				for descendant_filename in descendant_filenames:
					l+= 1
					print('Node ' + str(i) + ' out of ' + str(num_nodes) + ': ' + name + '; Descendant Node ' + str(k) + ' out of ' + str(num_descendants) + ': ' + descendant_name + '; Descendant File ' + str(l) + ' out of ' + str(descendant_numfiles))
					f = open(descendantfilenames_filename,'a')
					f.write(descendant_filename + '\n')
					f.close()

def read_bloomfiltersizes(bloomfiltersizes_filename): #The bloomfiltersizes file was created by running the KMC kmer counter on the files of descendant filenames created from the function write_descendantfiles (above)
	bloomfiltersizes = {}
	for line in open(bloomfiltersizes_filename):
		edited_name, size = line.strip().split(' ')
		bloomfiltersizes[edited_name] = int(size)
	return bloomfiltersizes

def construct_bloomfilter(bloomfiltersizes, num_nodes, nodesQueue, nodes_stop_event):
	ncbi = NCBITaxa()
	while not nodes_stop_event.is_set():
		i, node = nodesQueue.get()
		taxonid = int(node.name)
		name = ncbi.translate_to_names([taxonid])[0]
		print('Constructing bloom filter for node ' + str(i+1) + ' out of ' + str(num_nodes) + ': ' + name)
		edited_name = name.replace(' ', '_').replace('/', '_')
		#descendantfilenames_filename = edited_name + '.descendantfilenames'
		descendantfilenames_filename = edited_name + '.readfiles'
		bv_filename = edited_name + '.bv'
		if os.path.exists(bv_filename):
			print(bv_filename + ' already exists.')
			nodesQueue.task_done()
		else:
			size = bloomfiltersizes[edited_name] #the size was determined by counting unique kmers using kmc kmer counter
			size = int(round(size * 4.25,-3)) #multiply by 4.25; round to nearest 1000 (for a 13% false positive rate; bitvector length must be a multiple of 8)
			node.bf = BloomFilter(size, num_hashes = 3) #num_hashes = 3 (for a 13% positive rate)
			with open(descendantfilenames_filename) as f:
				num_files = sum(1 for _ in f)
			for j, descendant_filename in enumerate(open(descendantfilenames_filename)): #this for loop and the for loop below add all the kmers from all the descendant node files to the bloom filter
				#descendant_filename = descendant_filename.strip()
				descendant_filename = descendant_filename.strip() + '.dumps'
				print('Reading file ' + str(j+1) + ' out of ' + str(num_files) + ': ' + name)
				for line in open(descendant_filename):
					kmer = line.strip().split(' ')[0]
					node.bf.add(kmer)
			f = open(bv_filename, 'wb') #open bitvector file for writing
			node.bf.bv.tofile(f) #write bitvector to file
			f.close()
			print('Bloom filter for ' + name + ' constructed.')
			delattr(node, 'bf') #remove bloomfilter from memory
			nodesQueue.task_done()
	return

def construct_bloomfilters(tree, bloomfiltersizes):
	num_nodes = len(list(tree.get_descendants()))
	
	nodesQueue = queue.Queue()
	nodes_stop_event = threading.Event()
	
	for i, node in enumerate(tree.get_descendants()): #don't need to get a bloom filter for the root node 'Bacteria'
		nodesQueue.put((i,node))
	
	num_threads = 100
	threads = []
	for _ in range(num_threads):
		thread = threading.Thread(target = construct_bloomfilter, args = (bloomfiltersizes, num_nodes, nodesQueue, nodes_stop_event))
		thread.start()
		threads.append(thread)
	
	nodesQueue.join()
	nodes_stop_event.set()
	
	for t in threads:
		t.join()

def get_query_kmers(query, taxonid_to_readfilenames):
	ncbi = NCBITaxa()
	querykmers = []
	try: #if the query is a taxonid
		querytaxonid = int(query)
		queryname = ncbi.translate_to_names([querytaxonid])[0]
		print('Query name is ' + queryname)
		queryreadfilenames = taxonid_to_readfilenames[querytaxonid]
		for queryreadfilename in queryreadfilenames:
			querydumpsfilename = queryreadfilename + '.dumps'
			for line in open(querydumpsfilename):
				kmer = line.strip().split(' ')[0]
				querykmers.append(kmer)
	except: #if the query is a filename
		queryfilename = query
		print('Query filename is ' + queryfilename)
		for i, line in enumerate(open(queryfilename)):
			kmer = line.strip().split(' ')[0]
			querykmers.append(kmer)
			if i == 100000: #only count a given number of kmers, for testing purposes
				break
	return querykmers

def get_kmer_matches(next_node_kmers, childrenQueue, children_stop_event):
	ncbi = NCBITaxa()
	while not children_stop_event.is_set():
		if not childrenQueue.empty():
			child, current_kmers, threshold = childrenQueue.get()
			taxonid = int(child.name)
			name = ncbi.translate_to_names([taxonid])[0]
			edited_name = name.replace(' ', '_').replace('/', '_')
			bv_filename = edited_name + '.bv'
			print('Loading ' + name)
			child.bf = bf_from_bvfilename(bv_filename) #load bloom filter from bitvector
			print(name + ' loaded')
			kmer_matches = []
			for kmer in current_kmers: #figure out how many kmers of query match the current bloom filter
				if child.bf.contains(kmer):
					kmer_matches.append(kmer)
			delattr(child, 'bf') #remove bloom filter from memory
			if len(kmer_matches) > threshold: # if number of kmers that match exceeds threshold, add this child and matching kmers to the work queue
				next_node_kmers.append((child, kmer_matches))
			childrenQueue.task_done()
	return

def get_next_node_kmers(children, current_kmers, threshold):
	next_node_kmers = [] #a list of (node, kmers) tuples
	
	childrenQueue = queue.Queue()
	children_stop_event = threading.Event()
	
	for child in children:
		childrenQueue.put((child, current_kmers, threshold))
	
	num_threads = 100
	threads = []
	for _ in range(num_threads):
		thread = threading.Thread(target = get_kmer_matches, args = (next_node_kmers, childrenQueue, children_stop_event))
		thread.start()
		threads.append(thread)
	
	childrenQueue.join()
	children_stop_event.set()
	
	for t in threads:
		t.join()
	
	return next_node_kmers

def analyze_node(q, queried, queueLock, workQueue, threshold, num_kmers, responses, still_working, counterLock, queriedLock): #each thread targets this function to analyze its current node
	while not (still_working == 0 and workQueue.empty()): #while work is being done or the work queue is not empty
		queueLock.acquire()
		if not workQueue.empty():
			current_node, current_kmers = q.get()
			with counterLock:
				still_working += 1
			queueLock.release()
			ncbi = NCBITaxa()
			current_taxonid = int(current_node.name)
			current_name = ncbi.translate_to_names([current_taxonid])[0]
			print('Current node: ' + current_name)
			if current_name == 'Proteobacteria':
				with queueLock:
					for child in current_node.children:
						workQueue.put((child, current_kmers))
						print('adding')
					with counterLock:
						still_working -= 1
			elif current_node.is_leaf():
				responses[current_name] = len(current_kmers)/num_kmers
				print('Proportion of query kmers matching ' + current_name + ': ' + str(len(current_kmers)/num_kmers))
				with counterLock:
					still_working -= 1
			else:
				children = current_node.children
				node_kmers = get_next_node_kmers(children, current_kmers, threshold)
				with queriedLock:
					num_queried = queried.get()
					num_queried += len(children)
					queried.put(num_queried)
				if node_kmers == []:
					responses[current_name] = len(current_kmers)/num_kmers
					print('Proportion of query kmers matching ' + current_name + ': ' + str(len(current_kmers)/num_kmers))
					with counterLock:
						still_working -= 1
				else:
					with queueLock:
						for (child, kmer_matches) in node_kmers:
							workQueue.put((child, kmer_matches))
							print('adding')
						with counterLock:
							still_working -= 1
		else:
			with counterLock:
				still_working -= 1
			queueLock.release()
	return

#query the tree
def query_tree(tree, query, threshold, taxonid_to_readfilenames):
	
	#get query kmers and threshold
	querykmers = get_query_kmers(query, taxonid_to_readfilenames)
	num_kmers = len(querykmers)
	print('num_kmers', num_kmers)
	threshold = num_kmers * threshold
	
	#initialize the final variable responses (a dictionary with name as key and the percent of kmers matching the name as value)
	responses = {}
	
	#get the total number of nodes in the tree
	num_nodes = len(list(tree.traverse()))
	
	#create and initialize the queue and locks
	workQueue = queue.Queue() #this "work" queue contains all the nodes and corresponding matching kmers for each query path down the tree
	queueLock = threading.Lock()
	workQueue.put((tree,querykmers))
	
	queried = queue.Queue() #this queue counts how many nodes have been visited during the querying process
	queriedLock = threading.Lock()
	queried.put(0)
	
	#counterLock = threading.Lock()
	#still_working = 0 #this variable counts how many of the threads are currently doing work
	stop_event = threading.Event()
	
	#create new threads
	num_threads = 100
	threads = []
	for _ in range(num_threads):
		thread = threading.Thread(target = analyze_node, args = (workQueue,queried, queueLock, workQueue, threshold, num_kmers, responses, still_working, counterLock, queriedLock))
		thread.start()
		threads.append(thread)
	
	#wait for all threads to complete
	for t in threads:
		t.join()
	
	num_queried = queried.get()
	print(str(num_queried) + ' nodes queried out of ' + str(num_nodes))
	return sorted(responses.items(), key = operator.itemgetter(1), reverse = True)

if __name__=="__main__":
	
	command = sys.argv[1]
	
	num_taxonids = int(sys.argv[2])
	
	if command:
		taxonid_to_readfilenames = get_taxonid_to_readfilenames('name_ftpdirpaths')
		tree = get_tree(taxonid_to_readfilenames, num_taxonids)
		num_nodes = len(list(tree.traverse()))
		print('Tree has ' + str(num_nodes) + ' nodes')
	
	if command == "bloomfilters":
		#construct the bloomfilters (only necessary for the first time building the database)
		#actually, the end user never needs to perform this step, since they will download the bloom filters from the beginning
		#bloomfiltersizes = read_bloomfiltersizes('bloomfiltersizes')
		bloomfiltersizes = read_bloomfiltersizes('bloomfiltersizes2')
		construct_bloomfilters(tree, bloomfiltersizes)
	
	if command == "query":
		querytaxonid = sys.argv[3]
		threshold = float(sys.argv[4])
		matches = query_tree(tree, querytaxonid, threshold, taxonid_to_readfilenames)
		print(matches)
