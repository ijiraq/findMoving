vos_base_dir="vos:fraserw/LSST2020_June/psfOutputs"

function mcp {
    loop_counter=0
    while true
    do
        vcp $1 $2 && break || sleep 1 && loop_counter=$((loop_counter+1)) &&	[ "$loop_counter" -lt "3" ] || return
    done
}

for tarfile in $(vls ${vos_base_dir})
do
    mcp ${vos_base_dir}/${tarfile} ./
    visit=$(echo ${tarfile} | grep -Eo '[0-9]{6}')
    tar xf ${tarfile} && rm ${tarfile}
    ./fwhms.py $(find ./ -name "*${visit}*.psf.fits" -print)
    \rm -r HSC_June*
done
