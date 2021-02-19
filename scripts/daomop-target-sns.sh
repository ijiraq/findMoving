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
#  pointing   ccd   index         x         y   rate   angle    visit                                     stack                   ra                   dec   num
#     3071     2    1297    816.73    321.14    2.5     5.0   218194   STACK-0218194-002-00-+02.50-+05.00.fits    287.1405757175446   -20.816734590737564    3

basedir=$1 && shift
exptype=$1 && shift
pointing=$1 && shift
ccd=$1 && shift
index=$1 && shift
x=$1 && shift
y=$1 && shift
rate=$1 && shift
angle=$1 && shift
visit=$1 && shift
stack=$1 && shift
ra=$1 && shift
dec=$1 && shift
num=$1 && shift
if [ ${exptype} == "CORR" ]
then
    input="processCcdOutputs"
else
    input="diff"
fi
vosbase="vos:NewHorizons/S20A-OT04/"
output="sns" 
filter="HSC-R2"
section=750
chip=$(echo ${ccd} | awk '{printf("%03d",$1)}')
pointing=$(echo ${pointing} | awk '{printf("%05d",$1)}')
dbimages="vos:NewHorizons/dbimages/${pointing}/${chip}"

# Get from VOSpace the diff images to be stacked for this sources
\rm -r ${basedir}/rerun/diff/deepDiff
vcp  vos:NewHorizons/S20A-OT04/DIFFS/${pointing}/DIFF-${pointing}-${chip}_masked.tar ./
tar xvf DIFF-${pointing}-${chip}_masked.tar
rm DIFF-${pointing}-${chip}_masked.tar
for file in $(find ./HSC_All -name "*_masked.fits" -print)
do
    f=`echo ${file} | sed -e 's/_masked//'`
    mv ${file} ${f}
done


echo "Stacking ${section}X${section}  pixel box around ${ra} ${dec} at rate: ${rate} angle: ${angle}"
daomop-sns  ${basedir}\
    --swarp \
    --pointing ${pointing} \
    --rerun ${input}:${output} \
    --filter HSC-R2 \
    --ccd ${ccd} \
    --log-level INFO \
    --exptype ${exptype} \
    --group \
    --centre ${ra} ${dec} \
    --angle-min ${angle} \
    --angle-max ${angle} \
    --angle-step 1 \
    --rate-min ${rate} \
    --rate-step 1 \
    --rate-max ${rate} \
    --n-sub-stacks ${num} \
    --section-size ${section} \
    --clip 8

echo "Copying stacks to ${dbimages}"
for stack in ${basedir}/rerun/${output}/${exptype}/${pointing}/${filter}/*.fits
do
    expnum=$(basename ${stack} | awk -Fp '{printf("%s",$1)}')
    ccd=$(basename ${stack} | awk -Fp '{printf("%s",$2)}')
    ccd=$(echo ${ccd%.fits} | awk -Fp '{printf("ccd%02d", $1)}')
    vmkdir -p ${dbimages}/${expnum}/${ccd} || exit
    mcp ${stack} ${dbimages}/${expnum}/${ccd}/ || exit
    rm ${stack}
done

