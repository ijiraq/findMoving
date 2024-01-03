# setup aliases for daomop processing

next_cand() {
    d="$(pwd)"
    touch "${1}"
    rm MEASURE 
    cd ../ 
    [ -d "../${1}" ] || mkdir "../${1}"
    [ -f "../${1}" ] && echo "BAD DIR - ../${1}" && exit
    [ -d "../${1}/$(basename $d)" ] || mv "${d}" "../${1}/$(basename $d)"
    d="$(ls */MEASURE| head -n 1)"
    cd "${d%/MEASURE}"
    pwd
    ls
}

alias exam_cand='daomop-measure --dbimages . file track'
