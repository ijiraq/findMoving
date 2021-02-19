#!/bin/bash

function mcp {
    loop_counter=0
    while true
    do
	vcp -v $1 $2 && break || sleep 10 && loop_counter=$((loop_counter+1)) &&
 [ "$loop_counter" -lt "1000" ] || return 
    done
}

basedir=$1 && shift
pointing=$1 && shift
ccd=$1 && shift


vosbase="vos:NewHorizons/S20A-OT04/"
chip=$(echo ${ccd} | awk '{printf("%03d",$1)}')
pointing=$(echo ${pointing} | awk '{printf("%05d",$1)}')
dbimages="vos:NewHorizons/dbimages/${pointing}/${chip}"

# Get from VOSpace the diff images to be stacked for this sources
\rm -r ${basedir}/rerun/diff/deepDiff
vcp  vos:NewHorizons/S20A-OT04/DIFFS/${pointing}/DIFF-${pointing}-${chip}_masked.tar ./
tar xvf DIFF-${pointing}-${chip}_masked.tar
rm DIFF-${pointing}-${chip}_masked.tar


