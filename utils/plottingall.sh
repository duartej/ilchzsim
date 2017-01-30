#!/bin/bash 
# ------------------------------------------------------
# Plotting script (see `hzplots`).
# Note that the longitudinal impact parameter (-z option)
# is set to 40.0 mm, the maximum parallel momentum to scan
# is set to 40, and the parallel momentum cut is set to 
# square shape. 
# For each PID hypothesis, a folder is created.
# 
# To run it just :
# ./plottingall.sh [ OPTIONS ]
#
# Options:
#  -c --only-compare:  do not run fixed_pid but only collect
#                 the pkl files and run compare_pid
#
#
# jorge.duarte.campderros@cern.ch (2015-04-12)

# set default values
onlycompare=false

# parse command line arguments here
while [[ $# -gt 0 ]]
do
    key="$1"
    case $key in
	-c|--only-compare)
	    onlycompare=true
	    shift
	    ;;
	*)
	    ;;
    esac
    shift
done

n=0
#NCPU=`cat /proc/cpuinfo |grep processor|wc -l`

# collect the pkl files for the different PID modes
pklfiles=""

for i in noPID 010PID 005PID PID; 
do
    #if (("$NCPU" <= "$n"));
    #then
    #    wait;
    #    n=`ps aux|grep hzplots|wc -l`;
    #fi
    if [[ $onlycompare = false ]]
    then
	mkdir -p ${i} ;
	cd ${i};
	hzplots fixed_pid -s png -z 40.0 -p 40 --pLcut-type square ${i} ../processedhz_all.root #&
	cd -;
    fi
    pklfiles=$pklfiles" $i/d0cut_dict_$i.pkl"
    #n=$(($n+1))
done;
# Be sure everything finished
#wait;

# run hzplots in compare_pid mode
hzplots compare_pid -s png `echo $pklfiles`
