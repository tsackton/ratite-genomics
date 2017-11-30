module load BUSCO/3.0.2-fasrc01
export BUSCO_CONFIG_FILE=$(pwd)/config.ini
export AUGUSTUS_CONFIG_PATH="$HOME/sw/progs/augustus-3.0.3/config/"
for SEQ in $(ls *.fa);
do
	for LIN in aves_odb9 vertebrata_odb9;
	do
	echo "Running $SEQ with $LIN genes"
	run_BUSCO.py -i $SEQ -l $LIN -o ${SEQ}-${LIN} -m geno --cpu 16 -sp chicken &> ${SEQ}-${LIN}-BUSCO.log
	done
done
