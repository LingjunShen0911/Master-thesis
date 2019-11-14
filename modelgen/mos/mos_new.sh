#!/bin/bash

### check if saber was loaded ###
command -v saber >/dev/null 2>&1 || { echo >&2 "Please load Saber. Aborting.
Usage: module load saber."; exit 1; }
### check if matlab was loaded ###
command -v matlab >/dev/null 2>&1 || { echo >&2 "Please load MATLAB. Aborting.
Usage: module load matlab."; exit 1; }

### read config-file ###
### :: TODO :: ###
#n=0
#echo "read: config file:"
#for i in `cat ./MOS_config.dat | tr '.' '\n'` ; do
#   str=${str},${i}
#   let n=$n+1
#   var=`echo "var${n}"`
#   echo $var is ... ${i}
#done

### check args ###
if [[ $# -lt 10 ]]; then
    echo "Usage: <MODULNAME> <VGSMIN> <VGSMAX> <VGSSTEP> <VDSMIN> <VDSMAX> <VDSSTEP> <SIGMA> <matlab_only> <NameOfID>"
    echo "Example: ./mos_new.sh irff130 -4 8 0.48 -2 15 0.68 6 0 id"
    echo "Aborting."
    exit 1
fi

# set modul name 
MODULNAME="$1"
echo "device: $MODULNAME"


matlab_only="$9"
if [[ $matlab_only != 1 ]]; then
	# set simulation parameters					
	VGSMIN="$2" 
	echo "Vgs_min= $VGSMIN"
	VGSMAX="$3"
	echo "Vgs_max= $VGSMAX"
	VGSSTEP="$4"
	echo "Vgs_step= $VGSSTEP"

	VDSMIN="$5" 
	echo "Vds_min= $VDSMIN"
	VDSMAX="$6"
	echo "Vds_max= $VDSMAX"
	# set simulation timestep
	VDSSTEP="$7"
	echo "Vds_step= $VDSSTEP"

	# set number of breakpoints
	SIGMA="$8"
	echo "sigma= $SIGMA"

	# set name of drain current
	NameOfID="${10}"
	echo "name_of_ID= $NameOfID"

	rm -rf data/
	mkdir -p data/

	#replace mos in NETLIST.sin
	cp mos_template.sin data/netlist.sin
	sed "s/###MOS###/$MODULNAME/g" -i data/netlist.sin

	MACRO="data/netlist.ai_mcr"

	#create saber script
	echo "Guide:LoadDesign -design netlist.sin" >>"$MACRO"
	echo "Saber:Send {vary /v.vgs from $VGSMIN to $VGSMAX by $VGSSTEP }" >>"$MACRO"
	echo "Saber:Send {dc " >>"$MACRO"
	echo "dt (pins all,siglist /$MODULNAME.dut_mos/$NameOfID ,sweep /v.vds,swrange from $VDSMIN to $VDSMAX by $VDSSTEP}" >>"$MACRO"
	echo "Saber:Send {end}" >>"$MACRO"

	#save plotfile in .txt format 
	echo "set graph0 \$Graph(graph)" >> "$MACRO"
	echo "set pf1 [ScopeSigMgr:loadpffile netlist.dt.ai_pl]" >> "$MACRO"
	### catch Saber 2009.06 "$ScopeSigMgr(embedded): no such variable" error
	echo "catch { ScopeSigMgr:SaveAsText \$pf1 saber_output.txt }" >>"$MACRO"
	echo "exit" >>"$MACRO"

	# start the simulation in saber
	echo "start saber"
	cd data/
	aimsh -c "AimSession:SendShutdown" &> /dev/null
	saber -script ../"$MACRO"
	aimsh -c "AimSession:SendShutdown" &> /dev/null
	cd ..

	#convert .txt to .csv
	sed 's/\t/,/g'  data/saber_output.txt > data/saber_output.csv 

	#set number of breakpoints
	cp -f mos_new.m data/mos.m
	cp -f check_model.m data/check_model.m
	echo "element_name = '${MODULNAME}';" >> data/settings.m
	echo "reduce_factor = ${SIGMA};" >> data/settings.m
	echo "Vgs_min = ${VGSMIN};" >> data/settings.m
	echo "Vgs_max = ${VGSMAX};" >> data/settings.m
fi

#start MATLAB
echo "start MATLAB"
cd data/
matlab -nodisplay -nosplash -r 'try, mos, catch, exit(1), end; exit' > mos.log || { echo "ERROR: MATLAB failed!"; exit 1; }
cd ..

#replace model name in model.dat
if [ -f "data/model.xml" ]; then
	echo ""
	echo "DONE => 'model.xml' generated!"
	#show results
	#eog data/check_model.png &
	eog data/MOS_sampled.png &
	eog data/MOS_reduced.png &

	#echo file "MOS_warnings.log"
	while read LINE
	do
		echo $LINE
	done < data/MOS_warnings.log
else
	echo " "
	echo "ERROR: file ""model.dat"" for device $MODULNAME was NOT created! See also ""MOS_warnings.log""!"
fi


