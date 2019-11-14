#!/bin/bash

### check args ###
if [[ $# -lt 5 ]]; then
    echo "Usage: <MODULNAME> <VMIN> <VMAX> <VSTEP> <SIGMA>"
    echo "Aborting."
    exit 1
fi

### clean old data ###
rm -rf data/
mkdir -p data/

### set paramters ###
# set saber-modul name
MODULNAME="$1" #d1n4148
echo "bauelement $MODULNAME"
# set voltage range (Vmin, Vmax)
VMIN="$2" 
echo "Vd_min= $VMIN"
VMAX="$3"
echo "Vd_max= $VMAX"
# set voltage step size
VSTEP="$4"
echo "Step= $VSTEP"
# set number of breakpoints (sigma)
SIGMA="$5"
echo "sigma= $SIGMA"

###replace diode in NETLIST.sin ###
cp diode_template.sin data/netlist.sin
sed "s/^DIODE/$MODULNAME/" -i data/netlist.sin

MACRO="data/netlist.ai_mcr"

### create saber script ###
echo "Guide:LoadDesign -design netlist.sin" >>"$MACRO"
echo "Saber:Send {dc" >>"$MACRO"
echo "dt (pins all,siglist /$MODULNAME.dut_diode/id ,sweep /v_dc.ud,swrange from $VMIN to $VMAX by $VSTEP}" >>"$MACRO"

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
echo "saber done"

#convert .txt to .csv
#sed 's/\t/,/g'  data/saber_output.txt > data/saber_output.csv 

#set number of breakpoints
cp diode.m data/diode.m
#cp pwl_optimal_diode.m data/pwl_optimal_diode.m
cp optimalLinearApproximation.m data/optimalLinearApproximation.m
cp segmentErr.m data/segmentErr.m
echo "element_name = '${MODULNAME}';" >> data/settings.m
echo "sigma = ${SIGMA};" >> data/settings.m


#start MATLAB
#export MATLABPATH="${PWD}/data:${PWD}"
#echo "start MATLAB"
#cd data/
#matlab -nojvm -nosplash -r 'try, diode, catch, exit(1), end; exit' > diode.log || { echo "ERROR: MATLAB failed!"; exit 1; }
#matlab -nodisplay -nosplash -r 'try, diode, catch, exit(1), end; exit' > diode.log || { echo "ERROR: MATLAB failed!"; exit 1; }
#cd ..
#echo "one extra breakpoint was pasted"
#SIGMA_NEU= "$(($SIGMA+1))"
#echo "sigma= $SIGMA_NEU"
#echo "MATLAB done"

