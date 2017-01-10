import re, csv, sys, os
from ete3 import EvolTree

#test classes
ratites={'droNov', 'casCas', 'strCam', 'aptHaa', 'aptOwe', 'aptRow', 'rheAme', 'rhePen'}
bop={'aquChr', 'falPer', 'halLeu'}
wb={'anaPla', 'aptFor', 'balReg', 'chaVoc', 'egrGar', 'nipNip', 'pygAde', 'fulGla'}
vl={'calAnn', 'corBra', 'serCan', 'geoFor', 'melUnd', 'pseHum', 'taeGut', 'ficAlb'}
rand1={'colLiv', 'chaVoc', 'halLeu', 'taeGut', 'nipNip'}
rand2={'falPer', 'picPub', 'lepDis', 'melUnd', 'aquChr'}


#given name, test if ratite or not
def node_in_class(node,tree,testclass):
    #takes a node, an ete3 EvolTree, and a set defining a class of tips
    #returns true if i) the node is a leaf and is in the class set or
    #ii) the node is internal and all descendant leaves are in the class set
    #WARNING: assumes but does not check that nodenames are unique on the tree
    #WARNING: will return IndexError if nodename is not in tree
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

infile=sys.argv[1]

with open(infile) as f:
    lines=f.read().splitlines()


for line in lines:
    if line=="":
        continue
    else:
        fields=line.split("\t")
        hog=fields[1]
        tree=fields[0]
        try:
            t=EvolTree(fields[2])
        except:
            continue
        for node in t.traverse():
            #UGLY!
            isratite=node_in_class(node,t,ratites)
            isbop=node_in_class(node,t,bop)
            iswb=node_in_class(node,t,wb)
            isvl=node_in_class(node,t,vl)
            isr1=node_in_class(node,t,rand1)
            isr2=node_in_class(node,t,rand2)
            brstat=node.dist
            nname=node.name
            if nname=="":
                nname="-".join(node.get_leaf_names())
            try:
                pname=node.up.name
            except AttributeError:
                pname="root"
            if pname=="":
                pname="-".join(node.up.get_leaf_names())
            brname=pname + ":" + nname
            print(hog,tree,pname,nname,brname,brstat,isratite,isbop,iswb,isvl,isr1,isr2, sep=",", end="\n")
