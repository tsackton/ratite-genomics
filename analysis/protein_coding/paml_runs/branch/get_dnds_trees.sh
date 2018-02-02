#!/bin/bash

HOG=$1
SUBDIRNUM=`expr ${HOG##HOG} % 100`
printf -v SUBDIR "%04d" $SUBDIRNUM
PAMLDIR=$SUBDIR/$HOG	
FINALDIR=$PAMLDIR/$HOG.codeml.br.ctl.out
OUTFILE=$FINALDIR/br.out
echo $OUTFILE
DATE=$(date +%F)
#dn tree
DNTREE=$(grep -A1 "dN tree" $OUTFILE | perl -p -e "s/dN tree:\n/$HOG\t/g" | perl -p -e 's/--\n//g')
IFS=$'\n'; arrDN=($DNTREE); unset IFS
echo -e	"tree1\t${arrDN[0]}\n" >> dn_trees-$DATE.out
echo -e	"tree2\t${arrDN[1]}\n" >> dn_trees-$DATE.out
#ds tree
DSTREE=$(grep -A1 "dS tree" $OUTFILE | perl -p -e "s/dS tree:\n/$HOG\t/g" | perl -p -e 's/--\n//g')
IFS=$'\n'; arrDS=($DSTREE); unset IFS
echo -e "tree1\t${arrDS[0]}\n" >> ds_trees-$DATE.out
echo -e "tree2\t${arrDS[1]}\n" >> ds_trees-$DATE.out
#omega tree
WTREE=$(grep -A1 "w ratios as labels for TreeView:" $OUTFILE | perl -p -e "s/w ratios as labels for TreeView:\n/$HOG\t/g" | perl -p -e 's/--\n//g' | perl -p -e 's/\s+#/: /g')
IFS=$'\n'; arrW=($WTREE); unset IFS
echo -e "tree1\t${arrW[0]}\n" >> w_trees-$DATE.out
echo -e "tree2\t${arrW[1]}\n" >> w_trees-$DATE.out
