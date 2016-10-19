import os
import re
import sys
from collections import defaultdict

def getTime(path):
    # for a given path to an out file, checks to see if the time is given and returns that
    timepat = re.compile('^Time used:\s+(\d+:\d+:?\d*)\s*$', re.M)
    with open(path) as f:
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
                return "RUNNING"
        if totalSecs == 0:
            return "RUNNING"
        else:
            return totalSecs

with open("all_hogs") as hfile:
    hogs=[line.rstrip('\n') for line in hfile]

models=("br", "branch", "ancrec", "site", "branchsite", "branchsitenull")

status = defaultdict(dict)

for root, dirs, files in os.walk('.'):    
    for name in files:
        if name.endswith(".out"):
            filetime=getTime(os.path.join(root,name))
            analysis=name[:-4]
            try:
                hog=root.split(sep="/")[2]
            except IndexError:
                print("Couldn't get HOG id for", root, "skipping...", file=sys.stderr)
                next
            status[hog][analysis]=filetime

for hog in hogs:
    for model in models:
        if model in status.get(hog, {}):
            runtime = status[hog][model]
        else:
            runtime = "NOT_STARTED"
        print(hog, model, runtime, sep="\t")
