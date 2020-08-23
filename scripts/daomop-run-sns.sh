#!/bin/bash
source ${HOME}/.profile

wd=$(pwd)

# Update the software, if needed
cd ${HOME}/findMoving
git remote update
[ "$(git rev-parse master)" != "$(git rev-parse origin/master)" ] && git pull && python3 setup.py install --user 
export DAOMOP_VERSION=$(git rev-parse master)
echo ${DAOMOP_VERSION}
cd ${wd}

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
vos_uri=$1
shift
pointing=$1
shift
chip=$1
shift
ccd=$(echo ${chip} | awk ' { printf("%03d",$1) } ' )

dbimages="vos:NewHorizons/dbimages"
filter="HSC-R2"
exptype="deepDiff"
plantOutputs="plantOutputs"
processCcdOutputs="processCcdOutputs"
stack_rerun="sns_weighted"
diff_rerun="diff"

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

# pull the DIFF tar ball for this CCD from DIFFS/pointing directory on VOSpace.
sudo mkdir -p /mnt/${basedir} && sudo chown -R jkavelaars /mnt/${basedir} && ln -s /mnt/${basedir} ./ 
# mcp "${vos_uri}/DIFFS/${pointing}/DIFF-${pointing}-${ccd}_masked.*" ./  && tar xvf DIFF-${pointing}-${ccd}_masked.* 

# stack the data
echo "# uri:${vos_uri} pointing:${pointing} chip:${chip} amin:${angle_min} astep:${angle_step} amax:${angle_max} rmin:${rate_min} rstep:${rate_step} rmax:${rate_max}"
sync=$(which daomop-filesync.sh)

echo "Scheduling CRON to copy files to VOSPACE"
cron_file=/tmp/${RANDOM}
crontab -l > ${cron_file}
cmd="0 * * * * SHELL=/bin/bash && ${sync} $(pwd)/${basedir} ${dbimages} ${pointing} >> $HOME/sync.out 2>&1"
grep -q ${cmd} ${cron_file} ||  echo ${cmd} >> ${cron_file}
crontab ${cron_file}
crontab -l

daomop-sns ${basedir} \
	   --group \
	   --rerun ${diff_rerun}:${stack_rerun} \
	   --rate-min ${rate_min} \
	   --rate-max ${rate_max} \
	   --rate-step ${rate_step} \
	   --angle-min ${angle_min} \
	   --angle-max ${angle_max} \
	   --angle-step ${angle_step} \
	   --pointing ${pointing} \
	   --ccd ${chip} \
	   --mask \
	   --section-size 1200 \
	   --masked \
   	   --clip 8 \
	   --log-level INFO 

${sync} $(pwd)/${basedir} ${dbimages} ${pointing} 

