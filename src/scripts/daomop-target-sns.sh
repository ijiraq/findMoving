#!/bin/bash
# Some command defaults:
filter="HSC-R2"
loglevel="ERROR"

DEBUG="10"
INFO="20"
WARNING="30"
ERROR="40"

function logmsg() {
  msg_level=$(eval echo \$$1)
  log_level=$(eval echo \$"$loglevel")
  [ "${log_level}" -le "${msg_level}" ] && echo "${1}: ${2}"
  [ "${msg_level}" -ge "${ERROR}" ] && echo "EXIT CODE ${3}" && exit "${3}"
}

function show_help() {
cat << EOF
Usage: ${0##*/} -l loglevel base_dir exptype pointing chip index x y rate angle mag stack ra dec num

${0##*/} takes the inputs of where a source was found and creates a new set of stacks using a subset
aka cutout, of the pixel data centred around the RA/DEC provided at the rate/angle requested.
These cutouts can be used to measure the RA/DEC of the source in 'num' independent stacks of the
data taken on given night/field [aka pointing in LSST parlance]

The cutout stacked data will be written to: base_dir/rerun/TARGXXXXXX
Where XXXXXX is set by mktemp to temporarily store results which will
Ultimately end up in your current directory in a directory named pointing/chip/index
Expect that pointing/chip/index is unique.

options:
-l loglevel (DEBUG,INFO,ERROR) DEFAULT LEVEL: ${loglevel}
-h show_help

Where:
     base_dir: directory that contains the LSST Pipeline rerun directory with the images to stack.
     exptype : type of LSST Pipeline images to stack (normally DIFFEXP, sometimes CORR or MASKED) use UNKNOWN if you don't know.
     pointing: the LSST Pipeline pointing number for this data (e.g. 03148)
     chip    : the CCD to stack up
     index   : index number of source on that chip (used to create provisional name during measure step)
     x       : The x-pixel location of the source in the searched stack set (ignored)
     y       : The y-pixel location of the source in the searched stack set (ignored)
     rate    : rate of motion to stack at (-ve implies target moving west)
     angle   : Angle to stack at 0 degrees is West, 90 is North
     mag     : magnitude of the source, expected for artificial sources.
     stack   : name of the stack set where the srouce was found (ignored)
     ra      : RA to centre stack around
     dec     : DEC to centre stack around
     num     : How many groups of stacks to make, each group provides independent source measure
     name    : Provisional name of thing being stacked.
     epoch   : MJD of ra/dec position.
     dra     : RA ARC rate of motion of target at epoch (DeltaRA/cos(DEC))
     ddec    : DEC ARC Rate of motion of target at epoch
EOF
}

OPTIND=1
while getopts lhf: opt; do
  shift
  case "${opt}" in
    l) loglevel="${1}"
      shift
      ;;
    h) show_help
      exit 0
      ;;
    *) show_help >& 2
      exit 1
      ;;
  esac
done


# Need at least 3 arguments to proceed.
if [ $# -lt 3 ]; then
  show_help
  exit 1
fi

basedir=$1 && shift
exptypes=$1 && shift
input_filename=$1
if [ ! -f "${input_filename}" ]; then
  input_filename=$(mktemp "${TMPDIR:=/tmp}/DAOMOP_TARGET_XXXXXX")
  logmsg DEBUG "Creating temporary file ${input_filename} to hold input arguments"
  echo "${@}" >"${input_filename}"
fi

logmsg INFO "Reading inputs from ${input_filename}"

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
  provisional_name=$1 && shift
  epoch=$1 && shift
  dra=$1 && shift
  ddec=$1 && shift

  # Set the size of the cutout area to make in target mode.
  section=$(echo "${dra}" "${ddec}"| awk '{printf("%d",4*sqrt($1*$1+$2*$2)*3600/0.16)}')
  [ "${section}" -lt "100" ] && section=100
  [ "${section}" -gt "600" ] && section=600

  # put leading zeros in but remove them first, if they are already there.
  pointing=$(echo "${pointing}" | awk '{printf("%05d", $1)}')
  chip=$(echo "${ccd}" | awk ' {printf("%03d", $1)}')
  index=$(echo "${index}" | awk '{printf("%04d", $1)}')

  # echo "${basedir} ${exptype} ${pointing} ${chip} ${index} ${x} ${y} ${rate} ${angle} ${stack} ${ra} ${dec} ${num}"

  # check for MASKED or DIFFEXP in basedir/rerun/*/pointing and select input dir and exptype that is best suited.
  if [ "${exptypes}" == "UNKNOWN" ]; then
    exptypes="MASKED DIFFEXP CORR"
  fi
  for exptype in ${exptypes}; do
    logmsg DEBUG "Looking for ${exptype} files in ${basedir}/rerun with sub-directory ${pointing}"
    [ -d "${basedir}/rerun" ] || logmsg ERROR "RERUN directory ${basedir}/rerun does not exist" 1
    input=$(find "${basedir}/rerun" -path "*${pointing}*" -name "${exptype}-*-${chip}.f*" -print | head -n 1)
    logmsg DEBUG "Found this file: -->${input}<--"
    [ "${input}" ] || continue
    input="${input##*rerun/}"
    input="${input%%/*}"
    break
  done
  if [ ! "${input}" ]; then
    logmsg WARNING "FAILED TO GUESS INPUT RERUN going to next row"
    continue
  fi

  logmsg INFO "Using ${exptype} for ${pointing} in rerun directory ${input}"

  # Create locations to store the stamps... do this before we both making the stamps.
  stack_dir="${pointing}/${chip}/${index}"
  [ -d "${stack_dir}" ] || mkdir -p "${stack_dir}" || exit

  # Create temporary location that daomop-sns will store data to
  output=$(mktemp -d "${basedir}/rerun/TARGXXXXXX")
  output="${output##*rerun/}"
  lsst_dir="${basedir}/rerun/${output}/${exptype}/${pointing}/${filter}/"
  [ "$(ls -A ${lsst_dir} 2>/dev/null)" ] && logmsg ERROR "${lsst_dir} not empty" 2

  logmsg DEBUG "Using daomop-sns to create stack"
  logmsg INFO "Stacking ${section}X${section} pixel box around ${ra} ${dec} at rate: ${rate} angle: ${angle}"

  daomop-sns "${basedir}" \
    --swarp \
    --pointing "${pointing}" \
    --rerun "${input}":"${output}" \
    --filter "${filter}" \
    --ccd "${ccd}" \
    --log-level "${loglevel}" \
    --exptype "${exptype}" \
    --group \
    --centre "${ra}" "${dec}" "${epoch}" \
    --angle-min "${angle}" \
    --angle-max "${angle}" \
    --rate-min "${rate}" \
    --rate-max "${rate}" \
    --n-sub-stacks "${num}" \
    --section-size "${section}" \
    --clip 8

  logmsg INFO "Moving result from ${lsst_dir} to ${stack_dir}"
  # Move files out of LSST directory and into a stack directory in local FS and on VOS
  find "${lsst_dir}" -type f -name '*.fits' -exec mv {} "${stack_dir}" \;

done <"${input_filename}"
