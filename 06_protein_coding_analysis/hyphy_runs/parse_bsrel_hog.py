import re, csv, sys, os
from ete3 import EvolTree

#given name, test if ratite or not
def node_in_class(nodename,tree,testclass):
    #takes a node name, an ete3 EvolTree, and a set defining a class of tips
    #returns true if i) the node is a leaf and is in the class set or
    #ii) the node is internal and all descendant leaves are in the class set
    #WARNING: assumes but does not check that nodenames are unique on the tree
    #WARNING: will return IndexError if nodename is not in tree
    node=tree.search_nodes(name=nodename)[0]
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

#ugly, so ugly        
hog=sys.argv[1]    
treefile=sys.argv[2] 
resultsfile=sys.argv[3] 
test_to_use=sys.argv[4]

if not os.path.isfile(resultsfile):
    print(hog, 'NA', 'NA', 'NA', 'NA', sep="\t")
    quit()

if not os.path.isfile(treefile):
    print(hog, 'NA', 'NA', 'NA', 'NA', sep="\t")
    quit()

with open(treefile, 'r') as treefile:
    treestring=treefile.read().replace('\n', '')

treestring=re.sub(r"{\w+}", "", treestring)
t=EvolTree(treestring, format=1)
#ugly
ratites={'droNov', 'casCas', 'strCam', 'aptHaa', 'aptOwe', 'aptRow', 'rheAme', 'rhePen'}
vl={'calAnn', 'corBra', 'serCan', 'geoFor', 'melUnd', 'pseHum', 'taeGut', 'ficAlb'}
rand1={'colLiv', 'chaVoc', 'halLeu', 'taeGut', 'nipNip'}
rand2={'falPer', 'picPub', 'lepDis', 'melUnd', 'aquChr'}
if test_to_use == "ratites":
    testclade=ratites
elif test_to_use == "vl":
    testclade=vl
elif test_to_use == "rand1":
    testclade=rand1
elif test_to_use == "rand2":
    testclade=rand2
else:
   print("Select a valid test class!")
   quit()


#now read results into list
with open(resultsfile, 'r') as f:
    reader=csv.reader(f)
    res_list=list(reader)

#make output table, one line per hog, 
rs={'hog' : hog, 'target_sel' : 0, 'target_non' : 0, 'other_sel' : 0, 'other_non' : 0, 'target_sel_nom' : 0, 'other_sel_nom' : 0, 'nom_branches' : "", 'cor_branches' : ""}

for line in res_list:
    taxa=line[0]
    if t.search_nodes(name=taxa):
        pval=float(line[6])
        pvalholm=float(line[7])
        is_testclass=node_in_class(taxa, t, testclade)
        if pval < 0.05:
            if pvalholm < 0.05:
                if is_testclass:
                    rs['target_sel'] = rs['target_sel']+1
                else:
                    rs['other_sel'] = rs['other_sel']+1
                rs['cor_branches'] = rs['cor_branches'] + taxa + ":"
            if is_testclass:
                rs['target_sel_nom'] = rs['target_sel_nom']+1
            else:
                rs['other_sel_nom'] = rs['other_sel_nom']+1
            rs['nom_branches'] = rs['nom_branches'] + taxa + ":"
        else:
            if is_testclass:
                rs['target_non'] = rs['target_non']+1
            else:
                rs['other_non'] = rs['other_non']+1

print(rs['hog'], rs['target_sel'], rs['other_sel'], rs['target_sel_nom'], rs['other_sel_nom'], rs['target_non'], rs['other_non'], rs['cor_branches'][:-1], rs['nom_branches'][:-1], sep="\t")

