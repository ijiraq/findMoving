#!/bin/bash
source ${HOME}/.profile

function mcp {
    loop_counter=0
    while true
    do
	vcp -v $1 $2 && break || sleep 10 && loop_counter=$((loop_counter+1)) &&
 [ "$loop_counter" -lt "1000" ] || return 
    done
}

function gt { $(echo  ${1} ${2} | awk ' { if ( $1 > $2 ) { print("true") } else { print("false") }}' ); }

basedir=$1
shift
dbimages=$1
shift
pointing=$1
shift

stack_rerun="sns_weighted"
exptype="deepDiff"
filter="HSC-R2"

# stack the data

for filename in ${basedir}/rerun/${stack_rerun}/${exptype}/${pointing}/${filter}/*.fits
do
    output_filename=$(basename ${filename})
    storage=${dbimages}/$(echo ${output_filename} | awk -Fp ' { printf("%s/ccd%s",$1, $2) }')
    storage=${storage%.*}
    vmkdir -p ${storage}
    mcp ${filename} ${storage}/${output_filename}
done

