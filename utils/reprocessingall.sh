#!/bin/bash 
# ----------------------------------------------------------
# Create a thinned version of all the root files produced by
# the `ilchz` exec (using the simulateall.sh script) using
# all availables cores. All the output files produced are 
# finally merged into an unique one: `processedhz_all.root`
# 
# In the same folder where the output of the `ilchz` files 
# are, just :
# ./reprocessingall.sh
#
# jorge.duarte.campderros@cern.ch (2015-04-12)



FNAMEPROV="__PROVISIONAL__OUTPUT"
n=0
NCPU=`cat /proc/cpuinfo |grep processor|wc -l`
for i in `ls hz*_*.root`; 
do
    if (("$NCPU" <= "$n"));
    then
        wait;
        n=0;
    fi
    # Reprocessing
    processhzroot $i -o ${FNAMEPROV}_${RANDOM}.root &
    n=$(($n+1))
done;
# Be sure everything finished
wait;

# Merge all the files
hadd processedhz_all.root ${FNAMEPROV}_*.root && rm ${FNAMEPROV}_*.root;

