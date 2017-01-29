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
samples="uubar ddbar ssbar ccbar bbbar gg"
cmdfile=("1:onMode" "0:onMode" "2:onMode" "3:onMode" "4:onMode" "9:onMode")

n=0
for s in $samples; 
do
    mode=${cmdfile[$n]}
    sed -i.bak -e 's/onMode = [01]/onMode = 0/g' ilchz.cmnd 
    sed -i.bak "s/25:$mode = [01]/25:$mode = 1/g" ilchz.cmnd
    echo "Processing MODE:$mode"
    # Obtaining the Kaon-kaon background (assuming PID)
    ilchz ilchz.cmnd -f "kaons" -b -s 1 1000 -o hz${s}_PID_kaons_only.root &
    # Obtaining a 0.05 of mis-identification
    ilchz ilchz.cmnd -f "kaons_pions" -b -s 1 1000 -m 0.05 -o hz${s}_005PID_kaons_pions.root &
    # Obtaining a 0.05 of mis-identification
    ilchz ilchz.cmnd -f "kaons_pions" -b -s 1 1000 -m 0.10 -o hz${s}_010PID_kaons_pions.root &
    # Obtaining no PID (pions and kaons used indistinctly
    ilchz ilchz.cmnd -f "kaons_pions" -b -s 1 1000 -o hz${s}_noPID_kaons_pions.root &
    wait 
    n=$(($n+1))
done;

