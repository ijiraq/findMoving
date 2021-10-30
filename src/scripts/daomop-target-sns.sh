#!/bin/bash
# Some command defaults:
filter="HSC-R2"
loglevel="ERROR"

if [ $# -lt 3 ]; then
  echo "Usage: ${0} base_dir exptype pointing chip index x y rate angle mag stack ra dec num"
  echo ""
  echo "${0} takes the inputs of where a source was found and creates a new set of stacks using a subset,"
  echo " aka cutout, of the pixel data centred around the RA/DEC provided at the rate/angle requested."
  echo " These cutouts can be used to measure the RA/DEC of the source in 'num' independent stacks of the"
  echo " data taken on given night/field [aka pointing in LSST parlance]"
  echo ""
  echo " The cutout stacked data will be written to: base_dir/rerun/RANDOM "
  echo " Where RANDOM is created by mktemp to temporarily store results which will "
  echo " Ultimately endup in your current directory in a directory named pointing/chip/index "
  echo " Expect that pointing/chip/index is unique."
  echo ""
  echo " settings:"
  echo " loglevel: ${loglevel}"
  echo "   filter: ${filter}"
  echo ""
  echo "Where: "
  echo "       base_dir: directory that contains the LSST Pipeline rerun directory with the images to stack."
  echo "       exptype : type of LSST Pipeline images to stack (normally DIFFEXP, sometimes CORR or MASKED)."
  echo "       pointing: the LSST Pipeline pointing number for this data (e.g. 03148)"
  echo "       chip     : the CCD to stack up"
  echo "       index   : index number of source on that chip (used to create provisional name during measure step)"
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
exptypes=$1 && shift
input_filename=$1
if [ ! -f "${input_filename}" ]; then
  input_filename=$(mktemp "${TMPDIR:=/tmp}/XXXXXX_tract.txt")
  echo "${@}" >"${input_filename}"
  echo "Created temporary file ${input_filename} to hold input arguments"
fi

echo "Reading inputs from ${input_filename}"

while read -r line; do
  if [[ "${line}" =~ "pointing" ]]; then
    continue
  fi
  set ${line}
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

  # put leading zeros in but remove them first, if they are already there.
  pointing=$(echo "${pointing}" | awk '{printf("%05d", $1)}')
  chip=$(echo "${ccd}" | awk ' {printf("%03d", $1)}')
  index=$(echo "${index}" | awk '{printf("%04d", $1)}')

  echo "${pointing}.${chip}.${index}.${x}.${y}.${rate}.${angle}.${stack}.${ra}.${dec}.${num}"

  # check for MASKED or DIFFEXP in basedir/rerun/*/pointing and select input dir and exptype that is best suited.
  if [ "${exptypes}" == "UNKNOWN" ] ; then
    exptypes="MASKED DIFFEXP CORR"
  fi
  for exptype in ${exptypes}; do
    input=$(find "${basedir}/rerun" -path "*${pointing}*" -name "${exptype}-*-${chip}.f*" -print | head -n 1)
    [ "${input}" ] || continue
    input="${input##*rerun/}"
    input="${input%%/*}"
    break
  done
  if [ "${input}" ]; then
    echo "FAILED TO GUESS INPUT RERUN"
    exit 1 ;
  fi

  echo "Using ${exptype} for ${pointing} in rerun directory ${input}"


  section=200

  # Create locations to store the stamps... do this before we both making the stamps.
  stack_dir="${pointing}/${chip}/${index}"
  [ -d "${stack_dir}" ] || mkdir -p "${stack_dir}" || exit

  # Create temporary location that daomop-sns will store data to
  output=$(mktemp -d "${basedir}/rerun/XXXXXX")
  lsst_dir="${basedir}/rerun/${output}/${exptype}/${pointing}/${filter}/"
  [ "$(ls -A ${lsst_dir})" ] && echo "${lsst_dir} not empty, exiting" && exit

  echo "Using daomop-sns to create stack at ${lsst_dir}"
  echo "Stacking ${section}X${section} pixel box around ${ra} ${dec} at rate: ${rate} angle: ${angle}"
  echo "Then moving result from ${lsst_dir} to ${stack_dir}"

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

  # Move files out of LSST directory and into a stack directory in local FS and on VOS
  find "${lsst_dir}" -type f -name '*.fits' -exec mv {} "${stack_dir}" \;

done <"${input_filename}"
