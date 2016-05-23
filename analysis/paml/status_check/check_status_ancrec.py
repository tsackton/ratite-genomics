#!/usr/bin/env python

import os
import re
import sys
from collections import defaultdict

def uniqCount(filename):
    sofar = {}
    try:
        f = open(filename, 'r')
        for line in f:
           if line in sofar:
               next
           else:
               sofar[line] = 1
    
        lines = len(sofar.keys())
    except IOError:
        lines=0
    
    return lines

def getTime(path):
    # for a given path to an out file, checks to see if the time is given and returns that
    timepat = re.compile('^Time used:\s+(\d+:\d+:?\d*)\s*$', re.M)
    try:
        with open(os.path.join(path, "ancrec.out")) as f:
            output = f.read()
            times=timepat.findall(output)
            totalSecs = 0
        for tm in times:
            timeParts = [int(s) for s in tm.split(':')]
            if len(timeParts) == 3:
                totalSecs += (timeParts[0] * 60 + timeParts[1]) * 60 + timeParts[2]
            elif len(timeParts) == 2:
                totalSecs += (timeParts[0] * 60) + timeParts[1]
            elif len(timeParts) == 1:
                totalSecs = timeParts[0]
            else:
                return "FAILED"
        if totalSecs == 0:
            return "FAILED"
        else:
            return totalSecs
    except:
         return "FAILED"
         
def checkTree(path):
    #for a given path to an out file, get all the trees in the outfile
    treepat = re.compile('TREE\s+#\s+\d+')
    try:
        with open(os.path.join(path, "ancrec.out")) as f:
            output = f.read()
            trees=treepat.findall(output)
            return len(trees)
    except:
        return "NA"
        	
with open("all_hogs") as hfile:
    hogs=[line.rstrip('\n') for line in hfile]

for hog in hogs:
    toppath = '{:0>4}'.format(int(hog) % 100)
    # 0000/100/100.codeml.ancrec.ctl.out/
    fullpath = toppath + "/" + hog + "/" + hog + ".codeml.ancrec.ctl.out"
    filetime=getTime(fullpath)
    filetrees=checkTree(fullpath)
    treefile=fullpath + "/" + hog + ".final.nwk"
#    print(treefile)
    totaltrees=int(uniqCount(treefile))-1
    print(hog,filetime,filetrees,totaltrees,sep="\t")
