from make_tree import *
from construct_bloomfilters import *
from query_tree import *

tree = get_tree('name_ftpdirpaths')
tree_test = get_tree('name_ftpdirpaths', 10)

#construct the bloomfilters (only necessary for the first time building the database)
#actually, the end user never needs to perform this step, since they will download the bloom filters from the beginning
#construct_bloomfilters(tree)
#construct_bloomfilters(tree_test)
