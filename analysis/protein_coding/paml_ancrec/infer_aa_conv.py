from ete3 import EvolTree
import gzip, sys

#GLOBAL DEFS#
ratites={'droNov', 'casCas', 'strCam', 'aptHaa', 'aptOwe', 'aptRow', 'rheAme', 'rhePen'}

#FUNCTIONS#
def make_node_key(tree):
    #make lookup table of node id -> node
    nodes={}
    for node in t.traverse():
        nodes[node.node_id] = node
    return(nodes)

def parse_tree(treestring):
    t=EvolTree(treestring)
    return(t)

def node_in_class(nodename,tree,testclass):
    #takes a node name, an ete3 EvolTree, and a set defining a class of tips
    #returns true if i) the node is a leaf and is in the class set or
    #ii) the node is internal and all descendant leaves are in the class set
    #WARNING: assumes but does not check that nodenames are unique on the tree
    #WARNING: will return IndexError if nodename is not in tree
    try:
        node=tree.search_nodes(name=nodename)[0]
    except:
        node=nodename
        
    if node.is_leaf():
        if node.name in testclass:
            return True
        else:
            return False
    else:
        desc=set(node.get_leaf_names())
        check=desc - testclass
        if not check:
            return True
        else:
            return False

def get_mut_type(anc1,desc1,anc2,desc2):
    if anc1 == anc2:
        #start at same aa
        if desc1 == desc2:
            #end at same aa
            return("parallel")
        else:
            #end at different aa
            return("divergent")
    else:
        #start at different aa
        if desc1 == desc2:
            #end at same aa
            return("convergent")
        else:
            return("regular")

def test_ancestry(t,node1,node2):
    if node1 is node2:
        return("sister")
    elif node2 in node1.get_descendants():
        return("ancestral")
    elif node1 in node2.get_descendants():
        return("ancestral")
    else:
        return("clean")

def get_branch_type(t,node):
    if node.is_leaf():
        return("tip")
    else:
        return("internal")

def parse_muts(mutstring,t,nodes,testclass):
    muts=mutstring.split(",")
#   print(muts)
    pairs=[]
    for i in range(0,len(muts)):
#        print("i:", i, sep=" ")
        info1=muts[i].replace("..",":").replace("-",":").split(":")
        isratite1=node_in_class(nodes[int(info1[1])],t,testclass)
        brnames1="_".join(sorted(nodes[int(info1[1])].get_leaf_names()))
        for j in range(i+1,len(muts)):
#           print("j:", j, sep=" ")
            info2=muts[j].replace("..",":").replace("-",":").split(":")  
            brnames2="_".join(sorted(nodes[int(info2[1])].get_leaf_names())) 
            isratite2=node_in_class(nodes[int(info2[1])],t,testclass)
            branch_comp_label="-".join(sorted([nodes[int(info1[1])].get_topology_id(),nodes[int(info2[1])].get_topology_id()]))
            branch_names_label="-".join(sorted([brnames1,brnames2]))
            mut_type=get_mut_type(info1[2],info1[3],info2[2],info2[3])
            comp_type=int(isratite1)+int(isratite2)
#           print(int(info1[0]),int(info2[0]))
            try:
                ancestry1=nodes[int(info1[0])]
                ancestry2=nodes[int(info2[0])]
            except KeyError:
                print("ERROR at", info1, info2, sep="\t")
                break
            ancestry_type=test_ancestry(t,ancestry1,ancestry2)
            branch_comp_type="-".join(sorted([get_branch_type(t,nodes[int(info1[1])]),get_branch_type(t,nodes[int(info2[1])])]))
            pairs.append(",".join([str(branch_names_label),str(branch_comp_label),str(comp_type),str(ancestry_type),str(branch_comp_type),str(mut_type)]))
    return(pairs)


treefile=sys.argv[1]
mutfile=sys.argv[2]

#open tree file and parse
with gzip.open(treefile, 'rt') as f:
    treelines=f.read().splitlines()

trees={}
for line in treelines:
    fields=line.split()
    if fields[0] == "hog":
        continue
    else:
        try:
            trees[fields[0]][str(fields[2])]=fields[4]
        except KeyError:
            trees[fields[0]]={}
            trees[fields[0]][str(fields[2])]=fields[4]
            
#open mutfile and parse
results=[]
with gzip.open(mutfile, 'rt') as f:
    mutlines=f.read().splitlines()

curhog=''
curtree=''
ts=''
t=''
nodekey={}
for line in mutlines:
    fields=line.split()
    if (fields[0] == "hog"):
        continue
    else:
        prevhog=curhog
        prevtree=curtree
        curhog=fields[0]
        curtree=str(fields[1])
#       print(line, curhog, curtree, prevhog, prevtree)
        if (curhog != prevhog) or (curtree != prevtree):
            ts=trees[curhog][curtree]
            t=parse_tree(ts)
            nodekey=make_node_key(t)
        #now parse mut line
        for res in parse_muts(fields[3],t,nodekey,ratites):
            print(",".join([str(fields[0]),str(fields[1]),res]))

