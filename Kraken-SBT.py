from bitarray import bitarray
from collections import defaultdict
from ete3 import NCBITaxa
from multiprocessing import Pool, Value
import mmh3
import operator
import os
import queue
import sys
import threading

class BloomFilter:
	
	def __init__(self, size, num_hashes, bitvector = None):
		self.size = size
		self.num_hashes = num_hashes
		if bitvector is None:
			self.bv = bitarray(size)
			self.bv.setall(False)
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
def get_taxonid_to_readfilenames(name_ftpdirpaths_filename):
	ncbi = NCBITaxa()
	#search through all the data files and create a dictionary that maps taxonids to readfilenames
	taxonid_to_readfilenames = defaultdict(list) #a given taxonid may map to multiple readfiles, thus each value in the dictionary is a list
	for line in open(name_ftpdirpaths_filename):
		splits = line.strip().split(' ') #split on each of the spaces
		name = ' '.join(splits[:-1]) #concatenate everything before the last space
		ftpdirpath = splits[-1] #everything after the last space
		#readfilename = './Bacteria_Genomes/' + ftpdirpath.split('/')[-1] + '_genomic.fna.gz'
		readfilename = './Bacteria_Genomes/' + ftpdirpath.split('/')[-1] + '_genomic.fna'
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
				taxonid_to_readfilenames[taxonid].append(readfilename)
			if name == 'Mycobacterium intracellulare MOTT-64': #upon further inspection, this name in assembly_summary.txt is curiously absent from NCBI
				print(name+' is not in the NCBI database')
				sys.stdout.flush()
		else:
			taxonid = listtaxonid[0]
			taxonid_to_readfilenames[taxonid].append(readfilename)
	return taxonid_to_readfilenames

def bf_from_bvfilename(bv_filename):
	f = open(bv_filename, 'rb')
	bitvector = bitarray()
	bitvector.fromfile(f)
	f.close()
	bf = BloomFilter(bitvector.length(), 3, bitvector) #We use num_hashes = 3
	return bf

def get_tree(name_ftpdirpaths_filename,num_taxons = 0):
	ncbi = NCBITaxa()
	taxonid_to_readfilenames = get_taxonid_to_readfilenames(name_ftpdirpaths_filename)
	
	#get the desired number of unique taxonids, in order to create the phylogeny tree
	taxonids = taxonid_to_readfilenames.keys() #length is 6,740 as expected (6,741 minus the one not in NCBI)
	taxonids = list(set(taxonids)) #4,526 taxonids are unique
	if num_taxons != 0:
		taxonids = taxonids[:num_taxons] #smaller set of taxonids for tree construction and testing
	
	#return desired phylogeny tree
	return taxonid_to_readfilenames, ncbi.get_topology(taxonids) #5,360 total nodes for full dataset

def write_descendantfiles(tree):
	ncbi = NCBITaxa()
	i=0
	numnodes = len(list(tree.traverse()))
	for node in tree.traverse():
		i+=1
		taxonid = int(node.name)
		name = ncbi.translate_to_names([taxonid])[0]
		print('Node ' + str(i) + ' out of ' + str(numnodes) + ': ' + name)
		sys.stdout.flush()
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
				print('Node ' + str(i) + ' out of ' + str(numnodes) + ': ' + name + '; File ' + str(j) + ' out of ' + str(numfiles))
				sys.stdout.flush()
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
				print('Node ' + str(i) + ' out of ' + str(numnodes) + ': ' + name + '; Descendant Node ' + str(k) + ' out of ' + str(num_descendants) + ': ' + descendant_name)
				sys.stdout.flush()
				descendant_filenames = taxonid_to_readfilenames[descendant_taxonid]
				descendant_numfiles = len(descendant_filenames)
				l = 0
				for descendant_filename in descendant_filenames:
					l+= 1
					print('Node ' + str(i) + ' out of ' + str(numnodes) + ': ' + name + '; Descendant Node ' + str(k) + ' out of ' + str(num_descendants) + ': ' + descendant_name + '; Descendant File ' + str(l) + ' out of ' + str(descendant_numfiles))
					sys.stdout.flush()
					f = open(descendantfilenames_filename,'a')
					f.write(descendant_filename + '\n')
					f.close()

def read_bloomfiltersizes(bloomfiltersizes_filename):
	bloomfiltersizes = {}
	for line in open(bloomfiltersizes_filename):
		edited_name, size = line.strip().split(' ')
		bloomfiltersizes[edited_name] = int(size)
	return bloomfiltersizes
	
def construct_bloomfilters(tree, bloomfiltersizes):
	ncbi = NCBITaxa()
	i=0
	numnodes = len(list(tree.get_descendants()))
	for node in tree.get_descendants():
		i+=1
		taxonid = int(node.name)
		name = ncbi.translate_to_names([taxonid])[0]
		print('Node ' + str(i) + ' out of ' + str(numnodes) + ': ' + name)
		sys.stdout.flush()
		edited_name = name.replace(' ', '_').replace('/', '_')
		#descendantfilenames_filename = edited_name + '.descendantfilenames'
		descendantfilenames_filename = edited_name + '.readfiles'
		bv_filename = edited_name + '.bv'
		if False: #os.path.exists(bv_filename):
			print(bv_filename + ' already exists.')
			sys.stdout.flush()
		else:
			size = bloomfiltersizes[edited_name]
			size = int(round(size * 4.25,-3))
			node.bf = BloomFilter(size, num_hashes = 3)
			j = 0
			for line in open(descendantfilenames_filename):
				j += 1
				#descendant_filename = line.strip()
				descendant_filename = line.strip() + '.dumps'
				print('Reading file ' + str(j))
				for line2 in open(descendant_filename):
					kmer = line2.strip().split(' ')[0]
					node.bf.add(kmer)
			f = open(bv_filename, 'wb')
			node.bf.bv.tofile(f)
			f.close()
			delattr(node, 'bf')

def get_query_kmers(query, taxonid_to_readfilenames):
	ncbi = NCBITaxa()
	querykmers = []
	try:
		querytaxonid = int(query)
		queryname = ncbi.translate_to_names([querytaxonid])[0]
		print('Query name is ' + queryname)
		sys.stdout.flush()
		queryreadfilenames = taxonid_to_readfilenames[querytaxonid]
		for queryreadfilename in queryreadfilenames:
			querydumpsfilename = queryreadfilename + '.dumps'
			for line in open(querydumpsfilename):
				kmer = line.strip().split(' ')[0]
				querykmers.append(kmer)
	except:
		queryfilename = query
		print('Query filename is ' + queryfilename)
		sys.stdout.flush()
		i = 0
		for line in open(queryfilename):
			i += 1
			kmer = line.strip().split(' ')[0]
			querykmers.append(kmer)
			if i == 100000:
				break
	return querykmers

def get_kmer_matches(child, current_kmers, threshold):
	ncbi = NCBITaxa()
	taxonid = int(child.name)
	name = ncbi.translate_to_names([taxonid])[0]
	if name == 'Proteobacteria':
		return (child, current_kmers)
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
		return (child, kmer_matches)

def get_next_node_kmers(children, current_kmers, threshold):
	with Pool() as p:
		results = [item for item in p.starmap(get_kmer_matches, [(child, current_kmers, threshold) for child in children]) if item is not None]
	return results

def analyze_node(q, queried, queueLock, workQueue, threshold, num_kmers, responses, still_working, counterLock, queriedLock):
	while not (still_working == 0 and workQueue.empty()):
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
					with counterLock:
						still_working -= 1
			elif current_node.is_leaf():
				responses[current_name] = len(current_kmers)/num_kmers
				print('Proportion of query kmers matching ' + current_name + ': ' + str(len(current_kmers)/num_kmers))
				sys.stdout.flush()
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
					sys.stdout.flush()
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
	
	#initialize variables
	responses = {}
	num_threads = 5
	threads = []
	queueLock = threading.Lock()
	counterLock = threading.Lock()
	queriedLock = threading.Lock()
	queried = queue.Queue()
	queried.put(0)
	numnodes = len(list(tree.traverse()))
	
	#create and initialize the queue
	workQueue = queue.Queue()
	workQueue.put((tree,querykmers))
	still_working = 0
	
	#create new threads
	for _ in range(num_threads):
		thread = threading.Thread(target = analyze_node, args = (workQueue,queried, queueLock, workQueue, threshold, num_kmers, responses, still_working, counterLock, queriedLock))
		thread.start()
		threads.append(thread)
	
	#wait for all threads to complete
	for t in threads:
		t.join()
	
	num_queried = queried.get()
	print(str(num_queried) + ' nodes queried out of ' + str(numnodes))
	return sorted(responses.items(), key = operator.itemgetter(1), reverse = True)

if __name__=="__main__":
	
	#ncbi = NCBITaxa()
	
	command = sys.argv[1]
	
	num_taxons = int(sys.argv[2])
	
	if command:
		taxonid_to_readfilenames, tree = get_tree('name_ftpdirpaths', num_taxons)
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
