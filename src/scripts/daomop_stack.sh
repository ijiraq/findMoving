#!/bin/bash
# Given a candidate directory name stack that candidate based on MPC file present.
cand="${1}"

cd "${cand}"  || exit 1
find ./ -name "*.fits" -exec rm {} \;
for phase in TODO RESTACK REDO MEASURE STACKING track
do
   rm -f "${phase}"
done
[ "$(ls -A *.mpc 2>/dev/null)" ] || exit 1
[ "$(ls -A *.mpc | wc -l )" -gt "1" ] && exit 1

touch STACKING
daomop-predict-obs *.mpc track
daomop-target-sns.sh /arc/projects/NewHorizons/DATA/ UNKNOWN track && rm STACKING || touch SNS_FAIL

[ -f STACKING ] || touch MEASURE
