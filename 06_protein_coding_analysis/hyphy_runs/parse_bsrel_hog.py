import re, csv, sys, os
from ete3 import EvolTree

#given name, test if ratite or not
def node_in_class(nodename,tree):
    #takes a node name, an ete3 EvolTree, and a set defining a class of tips
    #returns true if i) the node is a leaf and is in the class set or
    #ii) the node is internal and all descendant leaves are in the class set
    #WARNING: assumes but does not check that nodenames are unique on the tree
    #WARNING: will return IndexError if nodename is not in tree
    node=tree.search_nodes(name=nodename)[0]
    if node.is_leaf():
        return node.name
    else:
        desc=sorted(set(node.get_leaf_names()))
        return "-".join(desc)
        
#ugly, so ugly        
hog=sys.argv[1]    
treefile=sys.argv[2] 
resultsfile=sys.argv[3] 

if not os.path.isfile(resultsfile):
    print(hog, 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', sep="\t")
    quit()

if not os.path.isfile(treefile):
    print(hog, 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', sep="\t")
    quit()

with open(treefile, 'r') as treefile:
    treestring=treefile.read().replace('\n', '')

treestring=re.sub(r"{\w+}", "", treestring)
t=EvolTree(treestring, format=1)

#now read results into list
with open(resultsfile, 'r') as f:
    reader=csv.reader(f)
    res_list=list(reader)

#make output table, one line per hog, 
rs={'hog' : hog, 'selected_nom' : 0, 'selected_holm' : 0, 'total_tests' : 0, 'nom_branches' : "", 'holm_branches' : "", 'tree' : treestring}

for line in res_list:
    taxa=line[0]
    if t.search_nodes(name=taxa):
        pval=float(line[6])
        pvalholm=float(line[7])
        node_id=trans_node(taxa, t)
        rs['total_tests'] = rs['total_tests'] + 1
        if pval < 0.05:
            if pvalholm < 0.05:
                rs['selected_holm'] = rs['selected_holm'] + 1
                rs['holm_branches'] = rs['holm_branches'] + node_id + ":"
            else:
                rs['selected_nom'] = rs['selected_nom'] + 1

print(rs['hog'], rs['selected_holm'], rs['selected_nom'], rs['total_tests'], rs['holm_branches'][:-1], rs['nom_branches'][:-1], rs['tree'], sep="\t")

