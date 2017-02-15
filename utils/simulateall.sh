#!/bin/bash 
# ------------------------------------------------------
# Script to simulate HZ-> qqbar using different particle
# identification detector (PID) hypothesis (a perfect PID
# a PID with 5% mis-id rate, a PID with 20% mis-id rate and 
# no PID at all). For each q=u,d,s,c,b,g the script runs 
# in parallel the four PID.
# 
# To run it just :
# ./simulateall.sh
#
# IMPORTANT: Do not change the script to run a different decay
#            mode in the same parallel loop, the needed input 
#            configuration file (ilchz.cmnd) is defining the
#            decay mode
#
# jorge.duarte.campderros@cern.ch (2015-04-12)



# variables
samples="ssbar ccbar bbbar gg"
cmdfile=("2:onMode" "3:onMode" "4:onMode" "9:onMode")

rmin=5
rmax=1000
n=0
for s in $samples; 
do
    mode=${cmdfile[$n]}
    sed -i.bak -e 's/onMode = [01]/onMode = 0/g' ilchz.cmnd 
    sed -i.bak "s/25:$mode = [01]/25:$mode = 1/g" ilchz.cmnd
    echo "Processing MODE:$mode"
    # perfect PID
    ilchz ilchz.cmnd -f "kaons" -b -s $rmin $rmax -e 1 0 1 -o hz${s}_1-0-1-PID_kaons_only.root &
    # realistic efficiencies for sqrt(2) separation
    ilchz ilchz.cmnd -f "kaons_pions" -b -s $rmin $rmax -e 0.9 0.44  0.75 -o hz${s}_0.9-0.44-0.75-PID_kaons_pions.root &
    ilchz ilchz.cmnd -f "kaons_pions" -b -s $rmin $rmax -e 0.8 0.28  0.75 -o hz${s}_0.8-0.28-0.75-PID_kaons_pions.root &
    ilchz ilchz.cmnd -f "kaons_pions" -b -s $rmin $rmax -e 0.8 0.28  0.8 -o hz${s}_0.8-0.28-0.8-PID_kaons_pions.root &
    ilchz ilchz.cmnd -f "kaons_pions" -b -s $rmin $rmax -e 0.7 0.19  0.75 -o hz${s}_0.7-0.19-0.75-PID_kaons_pions.root &
    #ilchz ilchz.cmnd -f "kaons_pions" -b -s $rmin $rmax -e 0.6 0.12  0.75 -o hz${s}_0.6-0.12-0.75-PID_kaons_pions.root &
    #ilchz ilchz.cmnd -f "kaons_pions" -b -s $rmin $rmax -e 0.4 0.05  0.75 -o hz${s}_0.4-0.05-0.75-PID_kaons_pions.root &
    #ilchz ilchz.cmnd -f "kaons_pions" -b -s $rmin $rmax -e 0.2 0.01  0.75 -o hz${s}_0.2-0.01-0.75-PID_kaons_pions.root &
    wait
    n=$(($n+1))
done;
