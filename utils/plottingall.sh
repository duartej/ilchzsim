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

for i in 1-1-1-PID 1-0-1-PID 0.8-0.3-0.75-PID 0.5-0.08-0.75-PID
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
	for j in 'CC' 'NC' 'NN' 'all'
	do
	    mkdir -p ${j}
	    cd ${j}
	    hzplots fixed_pid -s png -p 40 --pLcut-type square -R 5 -c ${j} ${i} ../../processedhz_all.root #&
	    cd ..
	done
	cd ..
    fi
    pklfilesall=$pklfilesall" ../$i/all/d0cut_dict_$i.pkl"
    pklfilesNN=$pklfilesNN" ../$i/NN/d0cut_dict_$i.pkl"
    pklfilesNC=$pklfilesNC" ../$i/NC/d0cut_dict_$i.pkl"
    pklfilesCC=$pklfilesCC" ../$i/CC/d0cut_dict_$i.pkl"
    #n=$(($n+1))
done;
# Be sure everything finished
#wait;

# run hzplots in compare_pid mode

for j in 'CC' 'NC' 'NN' 'all'
do
    mkdir -p ${j}
    cd ${j}
    pklfiles="pklfiles${j}"
    pkl=${!pklfiles}
    echo $pkl
    hzplots compare_pid -s png `echo ${pkl}`
    cd ..
done
