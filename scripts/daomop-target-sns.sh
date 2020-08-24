#!/bin/bash

function mcp {
    loop_counter=0
    while true
    do
	vcp -v $1 $2 && break || sleep 10 && loop_counter=$((loop_counter+1)) &&
 [ "$loop_counter" -lt "1000" ] || return 
    done
}

[ $# -eq 7 ] || echo "Usage: $0 {basedir} {pointing} {ccd} {ra} {dec} {rate} {angle}"
basedir=$1 && shift
pointing=$1 && shift
ccd=$1 && shift
ra=$1 && shift
dec=$1 && shift
rate=$1 && shift
angle=$1 && shift
vosbase="vos:NewHorizons/S20A-OT04/"
dbimages="vos:NewHorizons/dbimages/"
input="diff"
output="sns_cutout" 
source="DIFFS" 
filter="HSC-R2"
exptype="deepDiff"
chip=$(echo ${ccd} | awk '{printf("%03d",$1)}')

for tarfile in $(vls ${vosbase}/${source}/${pointing}/DIFF*${chip}*)
do 
    vopath=${vosbase}/${source}/${pointing}/${tarfile} 
    echo "copying ${vopath} to ${tarfile}"
    mcp ${vopath} ./${tarfile} || exit
    tar xvf ${tarfile} && rm ${tarfile} || exit
    echo "Stacking 100x100 pixel box around ${ra} ${dec} at rate: ${rate} angle: ${angle}"
    daomop-sns  ${basedir} \
		--pointing ${pointing} \
		--rerun ${input}:${output} \
		--filter HSC-R2 --ccd ${ccd} \
		--log-level INFO \
		--group \
		--centre ${ra} ${dec} \
		--angle-min ${angle} \
		--angle-max ${angle} \
		--angle-step 1 \
		--rate-min ${rate} \
		--rate-step 1 \
		--rate-max ${rate} \
		--section-size 100 \
		--clip 8

    echo "Copying stacks to ${dbimages}"
    for stack in ${basedir}/rerun/${input}/${exptype}/${pointing}/${filter}/*.fits
    do
	expnum=$(basename ${stack} | awk -Fp '{printf("%s",$1)}')
	ccd=$(basename ${stack} | awk -Fp '{printf("%s",$2)}')
	ccd=$(echo {ccd%.fits} | awk ' {printf("ccd$02d", $1)}')
	vmkdir -p ${dbimages}/${expnum}/${ccd} || exit
	mcp ${stack} ${dbimages}/${expnum}/${ccd}/ || exit
    done

    echo "Cleaning up the DIFFs"
    rm ${basedir}/rerun/${output}/${exptype}/${pointing}/${filter}/DIFF*${ccd}* || exit
done
