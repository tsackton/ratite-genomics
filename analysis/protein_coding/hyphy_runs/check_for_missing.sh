for HOG in $(cat all_hogs)
do
	SUBDIRNUM=`expr ${HOG##HOG} % 100`
	printf -v SUBDIR "%04d" $SUBDIRNUM
	RUNDIR=$SUBDIR/$HOG	
	OUT1=$RUNDIR/tree1/$HOG.aBSREL.OUT
	OUT2=$RUNDIR/tree2/$HOG.aBSREL.OUT
        TREE1=$RUNDIR/tree1/$HOG.tree1.nwk
	TREE2=$RUNDIR/tree2/$HOG.tree2.nwk

	if [[ ! -s $TREE2 ]]
	then
		cmp --silent $OUT1 $OUT2
		if [ $? -eq 0 ]
		then
			echo "$OUT2" >> tree2_to_remove
		fi
	fi
done
		 		
