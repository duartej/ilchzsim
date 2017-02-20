#!/bin/bash 
# ------------------------------------------------------
# Script to simulate HZ-> qqbar using different particle
# identification detector (PID) hypothesis (a perfect PID
# a PID with 5% mis-id rate, a PID with 20% mis-id rate and 
# no PID at all). For each q=u,d,s,c,b,g the script runs 
# in parallel the four PID.
# 
# To run it just :
# ./simulateall.sh [OPTIONS]
#
# OPTIONS:
#   -n  integer
#       number of events to be generated per channel [50000]
#   -r  float
#       minimal radius for Ks reconstruction [5]
#   -R  float
#       maximal radius for Ks reconstruction [1000]
#   -t  float
#       tracking efficiency for K, Ks and Pi, is multiplied with PID efficiency for K and Pi,
#       and reco efficiency for Ks [0.95]
#   -s  float
#       separation of the gaussians for K-pi discrimination [1.5]
#   -e  float
#       efficiency cut for K reconstruction, can be issued several times [0.7 0.8 0.9]
#   -k  float
#       efficiency for Ks reconstruction [0.9]
#
# IMPORTANT: Do not change the script to run a different decay
#            mode in the same parallel loop, the needed input 
#            configuration file (ilchz.cmnd) is defining the
#            decay mode
#
# jorge.duarte.campderros@cern.ch (2015-04-12)

# set default values

# number of events per channel
nevents=50000

# minimal and maximal radius for Ks reconstruction
rmin=5
rmax=1000

# tracking efficiency for the K, Ks and Pi
trackeff=0.95

# define the efficiency for Ks reconstruction and the separation of the gaussians in the kaon-pion
# distriction
efflist=(0.7 0.8 0.9)
kseff=0.9
separation=1.5


firsteff="1"
# parse command line arguments here
while getopts n:r:R:t:s:e:k: key
do
    case "${key}" in
	n)
	    nevents=${OPTARG}
	    ;;
	r)
	    rmin=${OPTARG}
	    ;;
	R)
	    rmax=${OPTARG}
	    ;;
	t)
	    trackeff=${OPTARG}
	    ;;
	s)
	    separation=${OPTARG}
	    ;;
	e)
	    if [[ "${firsteff}" == "1" ]]
	    then
		efflist=(${OPTARG})
		firsteff="0"
	    else
		efflist+=(${OPTARG})
	    fi
	    ;;
	k)
	    kseff=${OPTARG}
	    ;;
	*)
	    echo "unknown option $key"
	    exit
	    ;;
    esac
done

echo "    Nevents = ${nevents}"
echo "      r_min = ${rmin}"
echo "      r_max = ${rmax}"
echo "  trackeff. = ${trackeff}"
echo " separation = ${separation}"
echo "    efflist = ${efflist[@]}"
echo "      kseff = ${kseff}"

echo "" >> run.log
echo " =====================================" >> run.log
echo " parameters for the simulation:" >> run.log 
echo "    Nevents = ${nevents}" >> run.log
echo "      r_min = ${rmin}" >> run.log
echo "      r_max = ${rmax}" >> run.log
echo "  trackeff. = ${trackeff}" >> run.log
echo " separation = ${separation}" >> run.log
echo "    efflist = ${efflist[@]}" >> run.log
echo "      kseff = ${kseff}" >> run.log
echo " =====================================" >> run.log
echo "" >> run.log

# set the number of events in ilchz.cmd
sed -i.bak "s/Main:numberOfEvents.*! number of events to generate/Main:numberOfEvents = ${nevents} ! number of events to generate/g" ilchz.cmnd


# variables
samples="ssbar ccbar bbbar gg"
cmdfile=("2:onMode" "3:onMode" "4:onMode" "9:onMode")

n=0
for s in $samples; 
do
    mode=${cmdfile[$n]}
    sed -i.bak -e 's/onMode = [01]/onMode = 0/g' ilchz.cmnd 
    sed -i.bak "s/25:$mode = [01]/25:$mode = 1/g" ilchz.cmnd
    echo "Processing MODE:$mode"
    # perfect PID
    ilchz ilchz.cmnd -f "kaons" -b -s $rmin $rmax -t 0.95 -e 1 0 1 -o hz${s}_0.95-1-0-1-PID_kaons_only.root &

    for keff in "${efflist[@]}"
    do	
	pieff=`kaonpionseparation ${keff} ${separation}`
	echo "${keff}  ${pieff}"
	ilchz ilchz.cmnd -f "kaons_pions" -b -s $rmin $rmax -t $trackeff -e $keff $pieff $kseff -o hz${s}_${trackeff}-${keff}-${pieff}-${kseff}-PID_kaons_pions.root &
    done
    wait
    n=$(($n+1))
done;
