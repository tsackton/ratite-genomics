for RATITE in rheAme rhePen droNov casCas aptHaa aptRow aptOwe;
do
    for TREE in extended original;
    do
       sbatch launch_proteincoding_run.sh analyze_ratites.R $RATITE $TREE
       sbatch launch_proteincoding_run.sh analyze_flightless.R $RATITE $TREE
    done
done

for VL in ficAlb pseHum taeGut serCan geoFor corBra;
do
    for TREE in extended original;
    do
       sbatch launch_proteincoding_run.sh analyze_vocal_learners.R $VL $TREE
	done
done

for TREE in extended original;
do
	sbatch launch_proteincoding_run_array.sh analyze_random.R 3 $TREE
	sbatch launch_proteincoding_run_array.sh analyze_random.R 4 $TREE
done