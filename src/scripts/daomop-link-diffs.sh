#!/usr/bin/env bash
# Build a symlink tree of DIFFEXP images for daomop-sns (header-based stacking).
#
# Usage:
#   daomop-link-diffs.sh LISTFILE OUTDIR [--dbimages DIR] [--ccd FIRST-LAST | --ccd N,N,...] [--dry-run]
#
# LISTFILE lines (whitespace-separated; # comments and blanks skipped):
#   2773082/2773082p.fits 2022-08-01 AS1_July ...
# Uses: expnum (basename of col1 before /), night (col2), field (col3).
# Groups by (field, night); REFNUM = middle of sorted unique expnums (n//2).
# Links: OUTDIR/field/night/ccdXX/DIFFEXP-{EXPNUM}-{REFNUM}-{XX}.fits
#   ->   DBIMAGES/EXPNUM/ccdXX/DIFFEXP-{EXPNUM}-{XX}.fits
set -euo pipefail

usage() {
    cat <<EOF
Usage: $(basename "$0") LISTFILE OUTDIR [--dbimages DIR] [--ccd FIRST-LAST | --ccd N,N,...] [--dry-run]

  LISTFILE     Exposure list (field+night grouping from cols 1-3)
  OUTDIR       Root directory for symlink tree
  --dbimages   Source dbimages root (default: /arc/projects/classy/dbimages)
  --ccd        CCD range FIRST-LAST or comma list (default: 00-35)
  --dry-run    Print ln -s actions without creating links
EOF
    exit 1
}

DBIMAGES="/arc/projects/classy/dbimages"
CCD_SPEC="00-35"
DRY_RUN=0
LISTFILE=""
OUTDIR=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dbimages)
            [[ $# -ge 2 ]] || usage
            DBIMAGES="$2"
            shift 2
            ;;
        --ccd)
            [[ $# -ge 2 ]] || usage
            CCD_SPEC="$2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN=1
            shift
            ;;
        -h|--help)
            usage
            ;;
        -*)
            echo "Unknown option: $1" >&2
            usage
            ;;
        *)
            if [[ -z "${LISTFILE}" ]]; then
                LISTFILE="$1"
            elif [[ -z "${OUTDIR}" ]]; then
                OUTDIR="$1"
            else
                echo "Unexpected argument: $1" >&2
                usage
            fi
            shift
            ;;
    esac
done

[[ -n "${LISTFILE}" && -n "${OUTDIR}" ]] || usage
[[ -f "${LISTFILE}" ]] || { echo "LISTFILE not found: ${LISTFILE}" >&2; exit 1; }

# Expand --ccd into newline-separated zero-padded 2-digit CCD ids.
expand_ccd_spec() {
    local spec="$1"
    local part first last i padded
    local IFS=','
    # shellcheck disable=SC2086
    set -- ${spec}
    for part in "$@"; do
        part="${part//[[:space:]]/}"
        if [[ "${part}" == *-* ]]; then
            first="${part%%-*}"
            last="${part#*-}"
            if [[ ! "${first}" =~ ^[0-9]+$ || ! "${last}" =~ ^[0-9]+$ ]]; then
                echo "Invalid CCD range: ${part}" >&2
                exit 1
            fi
            first=$((10#${first}))
            last=$((10#${last}))
            if (( first > last )); then
                echo "Invalid CCD range (first > last): ${part}" >&2
                exit 1
            fi
            for (( i = first; i <= last; i++ )); do
                printf "%02d\n" "${i}"
            done
        elif [[ "${part}" =~ ^[0-9]+$ ]]; then
            printf "%02d\n" "$((10#${part}))"
        else
            echo "Invalid CCD spec: ${part}" >&2
            exit 1
        fi
    done | awk '!seen[$0]++'
}

CCD_LIST="$(expand_ccd_spec "${CCD_SPEC}")"
[[ -n "${CCD_LIST}" ]] || { echo "No CCDs selected from: ${CCD_SPEC}" >&2; exit 1; }

WORKDIR="$(mktemp -d "${TMPDIR:-/tmp}/daomop-link-diffs.XXXXXX")"
cleanup() { rm -rf "${WORKDIR}"; }
trap cleanup EXIT

# Parse list -> WORKDIR/groups as "field|night expnum" lines, then unique/sorted per group.
# Skip # comments and blank lines.
awk '
    {
        sub(/#.*/, "")
        if ($0 ~ /^[[:space:]]*$/) next
        if (NF < 3) {
            print "Warning: skipping malformed line (need >= 3 columns): " $0 > "/dev/stderr"
            next
        }
        path = $1
        night = $2
        field = $3
        # expnum = basename of path before /
        n = split(path, parts, "/")
        expnum = parts[1]
        if (expnum == "") {
            print "Warning: skipping line with empty expnum: " $0 > "/dev/stderr"
            next
        }
        print field "|" night, expnum
    }
' "${LISTFILE}" > "${WORKDIR}/raw.txt"

[[ -s "${WORKDIR}/raw.txt" ]] || {
    echo "No exposure groups found in ${LISTFILE}" >&2
    exit 1
}

# Unique (field|night, expnum), then sort by key then numeric expnum
sort -u -k1,1 -k2,2n "${WORKDIR}/raw.txt" > "${WORKDIR}/pairs.txt"

# Emit one line per group: KEY REFNUM EXPNUM1 EXPNUM2 ...
awk '
{
    key = $1
    en = $2
    if (key != prev && prev != "") {
        n = length(exps)
        # 0-based middle index: n//2
        ref = exps[int(n / 2)]
        printf "%s %s", prev, ref
        for (i = 0; i < n; i++) printf " %s", exps[i]
        printf "\n"
        delete exps
    }
    exps[length(exps)] = en
    prev = key
}
END {
    if (prev != "") {
        n = length(exps)
        ref = exps[int(n / 2)]
        printf "%s %s", prev, ref
        for (i = 0; i < n; i++) printf " %s", exps[i]
        printf "\n"
    }
}
' "${WORKDIR}/pairs.txt" > "${WORKDIR}/groups.txt"

[[ -s "${WORKDIR}/groups.txt" ]] || {
    echo "No exposure groups found after grouping" >&2
    exit 1
}

n_groups=0
n_links=0
n_missing=0

while read -r key REFNUM rest; do
    field="${key%%|*}"
    night="${key#*|}"
    # shellcheck disable=SC2086
    set -- ${rest}
    expnums=("$@")
    n=${#expnums[@]}

    n_groups=$((n_groups + 1))
    echo "Group ${field}/${night}: ${n} expnum(s), REFNUM=${REFNUM}"

    while IFS= read -r ccd; do
        [[ -z "${ccd}" ]] && continue
        dest_dir="${OUTDIR}/${field}/${night}/ccd${ccd}"
        if [[ "${DRY_RUN}" -eq 0 ]]; then
            mkdir -p "${dest_dir}"
        else
            echo "mkdir -p ${dest_dir}"
        fi

        for expnum in "${expnums[@]}"; do
            src="${DBIMAGES}/${expnum}/ccd${ccd}/DIFFEXP-${expnum}-${ccd}.fits"
            link_name="DIFFEXP-${expnum}-${REFNUM}-${ccd}.fits"
            dest="${dest_dir}/${link_name}"

            if [[ ! -f "${src}" ]]; then
                echo "Warning: missing source (skipped): ${src}" >&2
                n_missing=$((n_missing + 1))
                continue
            fi

            if command -v realpath >/dev/null 2>&1; then
                abs_src="$(realpath "${src}")"
            else
                abs_src="$(cd "$(dirname "${src}")" && pwd)/$(basename "${src}")"
            fi

            if [[ "${DRY_RUN}" -eq 1 ]]; then
                echo "ln -sfn ${abs_src} ${dest}"
            else
                ln -sfn "${abs_src}" "${dest}"
            fi
            n_links=$((n_links + 1))
        done
    done <<< "${CCD_LIST}"
done < "${WORKDIR}/groups.txt"

echo "----"
echo "Groups: ${n_groups}"
echo "Links created: ${n_links}"
echo "Missing sources: ${n_missing}"
if [[ "${DRY_RUN}" -eq 1 ]]; then
    echo "(dry-run: no files were created)"
fi
