#!/bin/bash
source ${HOME}/.profile

function mcp {
    loop_counter=0
    while true
    do
	vcp -v $1 $2 && break || (sleep 10 && loop_counter=$((loop_counter+1)) && [ "$loop_counter" -lt "1000" ] || return )
    done
}

function getcorr {
    listfile1=$1
    wc -l ${listfile1}
    ccd1=$(echo $2 | awk ' { printf("%d", $1) } ' )
    for lsstrun in "LSST2020" "LSST2020_June" "LSST2020_Aug"
    do
        base1="vos:NewHorizons/${lsstrun}/processCcdOutputs/"
	for expnum in $(cat ${listfile1})
	do
	    expnum=$(echo ${expnum} | awk ' { printf("%d", $1) } ')
	    tarfile="HSC_${expnum}_${ccd1}_output.tgz"
            mcp ${base1}/${tarfile} ./ || break
	done
	tar xvf ${tarfile} && rm ${tarfile} || return
    done
}

# cadc-get-cert --cert cadcproxy.pem || exit
[ $# -eq "3" ] || (echo "Usage: $0 HSC_DIR pointing ccd" ; exit )
basedir=$1
vos_uri="vos:NewHorizons/S20A-OT04/"
pointing=$2
chip=$3
expnum_list="${pointing}_expnum_list.txt"

filter="HSC-R2"
exptype="deepDiff"
plantOutputs="plantOutputs"
processCcdOutputs="processCcdOutputs"
stack_rerun="sns_weighted"
diff_rerun="diff"
ccd=$(echo ${chip} | awk ' { printf("%03d",$1) } ' )

# pull the DIFF tar ball for this CCD from DIFFS/pointing directory on VOSpace.
echo "Getting DIFFS and making MASKED version"
# get the psf file for this run
mcp "${vos_uri}/psf_fwhm.txt" ./ || exit

# pull the DIFF tar ball if _masked didn't exist for this CCD from DIFFS/pointing directory on VOSpace.
mcp "${vos_uri}/DIFFS/${pointing}/DIFF-${ccd}.*" ./  && tar xvf DIFF-${ccd}.* | awk -F- '{print($4)}' > ${expnum_list}  || exit

# pull the CORR tar ball for this CCD from DIFFS/pointing directory on VOSpace.
# mcp "${vos_uri}/CORR/${pointing}/CORR-${ccd}.*" ./ && tar xvf CORR-${ccd}.* && rm CORR-${ccd}.*  ||
getcorr ${expnum_list} ${ccd} || exit

# pull the MASK tar ball for this CCD from DIFFS/pointing directory on VOSpace to only redo missing ones.
mcp "${vos_uri}/DIFFS/${pointing}/MASKED-${ccd}.*" ./ && tar xvf MASKED-${ccd}.* && rm MASKED-${ccd}.*  || exit

# Build the MASKED versions of the files.
for fullpath in ${basedir}/rerun/${diff_rerun}/${exptype}/${pointing}/${filter}/DIFFEXP-*-${ccd}.fits
do
    echo ${fullpath}
    masked="${fullpath/DIFFEXP/MASKED}"
    echo ${masked}
    [ -f ${masked} ] && continue
    diff_filename=$(basename ${fullpath})
    visit=$(echo ${diff_filename} | grep -Eo '[0-9]{7}')
    echo ${visit}
    #visit=${visit##0}
    fwhm=$(grep "${visit#0} ${ccd}" psf_fwhm.txt | awk ' { print $3 } ')
    [ "Z"${fwhm} == "Z" ] && fwhm=5.0
    daomop-intelligentMasker ${basedir} \
			     --rerun ${processCcdOutputs}:${diff_rerun} \
			     --pointing ${pointing} \
			     --ccd ${chip} \
			     --visit ${visit} \
			     --clip 12 \
			     --pad 3 \
			     --psf-fwhm ${fwhm}  \
			     --log-level INFO
done

echo "Putting MASKED versions to VOSpace"
tarfile=MASKED-${ccd}.tar
masked_dir=${basedir}/rerun/${diff_rerun}/${exptype}/${pointing}/${filter}
tar czf ${tarfile} ${masked_dir}/MASKED*${ccd}.fits  || exit
echo "Copying ${tarfile} to VOSpace"
vmkdir -p ${vos_uri}/DIFFS/${pointing}
mcp ${tarfile} ${vos_uri}/DIFFS/${pointing}/${tarfile} ||  exit
echo "Cleaing up by deleting ${tarfile}"
rm ${tarfile} 
