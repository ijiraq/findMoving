#!/bin/bash

# Some command defaults:
output="calib"
filter="HSC-R2"
loglevel="ERROR"

# Handy function to repeatedly attempt to put data to VOSpace.
function mcp() {
    loop_counter=0
    while true; do
        vcp -v "$1" "$2" && break || sleep 10 && loop_counter=$((loop_counter + 1)) &&
          [ "${loop_counter}" -lt "10" ] || return
    done
}


if [ $# -lt 3 ]; then
    echo "Usage: ${0} base_dir exptype pointing ccd index x y rate angle mag stack ra dec num"
    echo ""
    echo "${0} takes the inputs of where a source was found and creates a new set of stacks using a subset,"
    echo " aka cutout, of the pixel data centred around the RA/DEC provided at the rate/angle requested."
    echo " These cutouts can be used to measure the RA/DEC of the source in 'num' independent stacks of the"
    echo " data taken on given night/field [aka pointing in LSST parlance]"
    echo ""
    echo " The cutout stacked data will be written to: base_dir/rerun/${output} "
    echo " settings:"
    echo " loglevel: ${loglevel}"
    echo "   output: ${output}"
    echo "   filter: ${filter}"
    echo ""
    echo "Where: "
    echo "       base_dir: directory that contains the LSST Pipeline rerun directory with the images to stack."
    echo "       exptype : type of LSST Pipeline images to stack (normally DIFFEXP, sometimes CORR or MASKED)."
    echo "       pointing: the LSST Pipeline pointing number for this data (e.g. 03148)"
    echo "       ccd     : the CCD to stack up"
    echo "       index   : index number of source on that ccd (used to create provisional name during measure step)"
    echo "       x       : The x-pixel location of the source in the searched stack set (ignored)"
    echo "       y       : The y-pixel location of the source in the searched stack set (ignored)"
    echo "       rate    : rate of motion to stack at (-ve implies target moving west)"
    echo "       angle   : Angle to stack at 0 degrees is West, 90 is North"
    echo "       mag     : magnitude of the source, expected for artificial sources."
    echo "       stack   : name of the stack set where the srouce was found (ignored)"
    echo "       ra      : RA to centre stack around"
    echo "       dec     : DEC to centre stack around"
    echo "       num     : How many groups of stacks to make, each group provides independent source measure"
    exit 0
fi

basedir=$1 && shift
exptype=$1 && shift
input_filename=$1
if [ ! -f "${input_filename}" ]; then
    input_filename="${TMPDIR:-/tmp}/track_file.txt"
    echo "${@}" >"${input_filename}"
    echo "Created temporary file ${input_filename} to hold input arguments "
fi

echo "Reading inputs from ${input_filename}"

while read -r line; do
    set ${line}
    if [ "$1" == "#" ]; then
        continue
    fi
    pointing=$1 && shift
    ccd=$1 && shift
    index=$1 && shift
    x=$1 && shift
    y=$1 && shift
    rate=$1 && shift
    angle=$1 && shift
    mag=$1 && shift
    stack=$1 && shift
    ra=$1 && shift
    dec=$1 && shift
    num=$1 && shift
    echo "${pointing}.${ccd}.${index}.${x}.${y}.${rate}.${angle}.${stack}.${ra}.${dec}.${num}"
    if [ "${exptype}" == "CORR" ]; then
        input="processCcdOutputs"
    else
        input="diff"
    fi

    section=500
    # put leading zeros in but remove them first, if they are already there.
    chip=$(printf %03d "${ccd##0}")
    pointing=$(printf %05d "${pointing##0}")
    index=$(printf %04d "${index##0}")

    # Create locations to store the stamps... do this before we both making the stamps.
    dbimages="vos:NewHorizons/dbimages/"
    stack_dir="${pointing}/${chip}/${index}"
    [ -d "${stack_dir}" ] || mkdir -p "${stack_dir}" || exit
    vmkdir -p "${dbimages}/${stack_dir}" || exit
    echo "Storing stacked stamps in ${stack_dir} and ${dbimages}/${stack_dir}"

    echo "Stacking ${section}X${section}  pixel box around ${ra} ${dec} at rate: ${rate} angle: ${angle}"
    daomop-sns "${basedir}" \
        --swarp \
        --pointing "${pointing}" \
        --rerun "${input}":"${output}" \
        --filter "${filter}" \
        --ccd "${ccd}" \
        --log-level "${loglevel}" \
        --exptype "${exptype}" \
        --group \
        --centre "${ra}" "${dec}" \
        --angle-min "${angle}" \
        --angle-max "${angle}" \
        --rate-min "${rate}" \
        --rate-max "${rate}" \
        --n-sub-stacks "${num}" \
        --section-size "${section}" \
        --clip 8

    lsst_dir="${basedir}/rerun/${output}/${exptype}/${pointing}/${filter}/"
    # Move files out of LSST directory and into a stack directory in local FS and on VOS
    for stack in "${lsst_dir}"*.fits; do
        mv "${stack}" "${stack_dir}"
        stack="${stack_dir}/"$(basename "${stack}")
        mcp "${stack}" "${dbimages}/${stack_dir}" || exit
    done

done <"${input_filename}"
