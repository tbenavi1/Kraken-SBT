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

if __name__=="__main__":
  from make_tree import *
  tree = get_tree('name_ftpdirpaths')
  tree_test = get_tree('name_ftpdirpaths',10)
  construct_bloomfilters
