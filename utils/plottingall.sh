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
# ./plottingall.sh
#
#
# jorge.duarte.campderros@cern.ch (2015-04-12)


n=0
#NCPU=`cat /proc/cpuinfo |grep processor|wc -l`

# collect the pkl files for the different PID modes
pklfiles=""

for i in noPID 020PID 005PID PID; 
do
    #if (("$NCPU" <= "$n"));
    #then
    #    wait;
    #    n=`ps aux|grep hzplots|wc -l`;
    #fi
    mkdir -p ${i} ;
    cd ${i};
    hzplots fixed_pid -s png -z 40.0 -p 40 --pLcut-type square ${i} ../processedhz_all.root #&
    pklfiles=$pklfiles" $i/d0cut_dict_$i.pkl"
    cd -;
    #n=$(($n+1))
done;
# Be sure everything finished
#wait;

# run hzplots in compare_pid mode
hzplots compare_pid -s png `echo $pklfiles`
