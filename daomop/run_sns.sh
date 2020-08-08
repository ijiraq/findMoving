#!/bin/bash

source ${HOME}/.profile
wd=$(pwd)
git clone https://github.com/ijiraq/findMoving
cd findMoving
python3 setup.py install --user
cd $wd

function mcp {
    loop_counter=0
    while true
    do
	vcp -v $1 $2 && break || sleep 10 && loop_counter=$((loop_counter+1)) && [ "$loop_counter" -lt "1000" ] || return 
    done
}

# cadc-get-cert --cert cadcproxy.pem || exit
vos_uri=$1
shift
pointing=$1
shift
chip=$1
shift

basedir=$(pwd)"/HSC_June19-lsst"
filter="HSC-R2"
exptype="deepDiff"
stack_rerun="sns_weighted"
diff_rerun="diff"
ccd=$(echo ${chip} | awk ' { printf("%03d",$1) } ' )

if [ $# -eq 2 ] 
then
    angle_min=$1 && shift
    rate_min=$1 && shift
    angle_max=${angle_min}
    rate_max=${rate_min}
    angle_step=1
    rate_step=1
else
    angle_min=-10 && [ "$#" -gt "0" ] && angle_min=$1 && shift
    angle_step=2.5 && [ "$#" -gt "0" ] && angle_step=$1 && shift
    angle_max=10 && [ "$#" -gt "0" ] && angle_max=$1 && shift
    rate_min=0.5 && [ "$#" -gt "0" ] && rate_min=$1 && shift
    rate_step=0.5 && [ "$#" -gt "0" ] && rate_step=$1 && shift
    rate_max=3.5 && [ "$#" -gt "0" ] && rate_max=$1 && shift
fi

storage=${vos_uri}/STACKS_V2/${pointing}/${ccd}
echo "# uri:${vos_uri} pointing:${pointing} chip:${chip} amax:${angle_min} astep:${angle_step} amin:${angle_min} rmin:${rate_min} rstep:${rate_step} rmax:${rate_max}"
vmkdir -p ${storage}
for rate in $(seq ${rate_min} ${rate_step} ${rate_max});
do
    echo "${rate}"
    for angle in $(seq ${angle_min} ${angle_step} ${angle_max})
    do
	echo "${angle}"
	echo "${pattern}"
	pattern=$(echo ${chip} ${rate} ${angle} | awk ' { printf("-%03d-0[0-2]-%+06.2f-%+06.2f.fits.fz", $1, $2, $3 ) } ')
	[ $(vls ${storage} | grep -- ${pattern} | wc -l) -eq 3 ] && echo "${pattern} already done " && continue 
	mcp ${vos_uri}/DIFFS/${pointing}/DIFF-${ccd}.tar ./ 
	tar xvf DIFF-${ccd}.tar || exit
	daomop-sns ${basedir} \
		   --rerun ${diff_rerun}:${stack_rerun} \
		   --rate-min ${rate} \
		   --rate-max ${rate} \
		   --rate-step ${rate_step} \
		   --angle-min ${angle} \
		   --angle-max ${angle} \
		   --angle-step ${angle_step} \
		   --pointing ${pointing} \
		   --ccd ${chip} \
		   --mask \
		   --section-size 600 \
   		   --clip 16 \
	       --log-level INFO  || exit
	mcp "${basedir}/rerun/${stack_rerun}/${exptype}/${pointing}/${filter}/STACK*.fits.fz" ${storage}/ || exit
	rm ${basedir}/rerun/${stack_rerun}/${exptype}/${pointing}/${filter}/STACK*.fits.fz || exit
   done
done
