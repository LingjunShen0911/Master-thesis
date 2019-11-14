#!/bin/bash

#########################################################
### CAUTION: Existing results will be overwritten !!! ###
#########################################################

#make results directory
rm -rf data
mkdir -p results

MODULNAME="irff130"

for SIGMA in $(seq 1 10); do

	#SIGMA=5

	UGSMIN=0
	UGSMAX=20
	printf -v UGSSTEP '%0.5f' "$(echo "($UGSMAX - $UGSMIN) / $SIGMA" | bc -l)"
	UDSMIN=0
	UDSMAX=20
	printf -v UDSSTEP '%0.5f' "$(echo "($UDSMAX - $UDSMIN) / $SIGMA" | bc -l)"

	#create results
	for i in $(seq 1 10); do

		printf -v REDUCE_FACTOR '%0.2f' "$(echo "${i} / 10" | bc -l)"

		echo " "
		echo "==> creating MOSFET model $MODULNAME with REDUCE_FACTOR=$REDUCE_FACTOR , SIGMA=$SIGMA"

		#make results directory for current model
		RESULT_DIRECTORY="$MODULNAME"_"$SIGMA"_"$REDUCE_FACTOR"
		mkdir -p results/"$RESULT_DIRECTORY"/

		#run mos.sh script
		echo "./mos.sh $MODULNAME $UGSMIN $UGSMAX $UGSSTEP $UDSMIN $UDSMAX $UDSSTEP $REDUCE_FACTOR 0"
		./mos.sh "$MODULNAME" $UGSMIN $UGSMAX $UGSSTEP $UDSMIN $UDSMAX $UDSSTEP $REDUCE_FACTOR 0 > console.log
		mv -f console.log results/"$RESULT_DIRECTORY"/console.log

		#copy results
		mv -f data/*.* results/"$RESULT_DIRECTORY"/
		echo " "
	done
done

echo "DONE!"
echo ""

