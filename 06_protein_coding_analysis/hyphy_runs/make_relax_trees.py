from ete3 import EvolTree
import sys

treepath=sys.argv[1]
treeout=sys.argv[2]

t = EvolTree(treepath)
ratites = {'aptHaa', 'aptRow', 'aptOwe', 'strCam', 'droNov', 'casCas', 'rheAme', 'rhePen'}
#annotate leaves
for leaf in t.traverse():
    if leaf.is_leaf():
        if leaf.name in ratites:
            leaf.add_features(mark="{RatiteLeaf}")
    else:
        #internal node, get all leaf names and make sure all are in ratites
        desc=set(leaf.get_leaf_names())
        checkDesc=desc - ratites
        if not checkDesc:
            leaf.add_features(mark="{RatiteInternal}")

t.write(outfile=treeout)
