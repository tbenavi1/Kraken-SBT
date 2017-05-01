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
