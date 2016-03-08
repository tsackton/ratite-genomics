#run dless

EXPLEN=45
TARCOV=0.30
VER=$1

#set up
for INPUT in $(ls /n/regal/edwards_lab/ratites/wga/phast/final_mafs/split_all/*.ss)
do
	BASE=${INPUT##*/}
	SAMP=${BASE%%.ss}
	phastCons --expected-length=$EXPLEN --target-coverage=$TARCOV --rho 0.3168 --msa-format SS $INPUT ../neutMods/${VER}_final.named.mod > $VER.dless.gff &
done
cd ..
