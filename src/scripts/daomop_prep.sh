#!/bin/bash
# Prepare a candidate directory for to run findMoving candidate vetting


cand="${1}"
echo "Setting up directory for ${cand}"
basedir="/arc/projects/NewHorizons/DATA/"

# Don't run prep if there is already vetting marker files.
[ -f "$1/STACKING" ] && exit
[ -f "$1/MEASURE" ] && exit
[ -f "$1/TODO" ] && exit

cd "${cand}" || exit 1
# Clean up links that might already be here.
rm -f *.reg
rm -f *.csv
ln -s ${basedir}/PlantList/ds9_reg_files/*.reg ./

ln -s ${basedir}/skycat.csv ./

rm -f track
touch TODO
