#!/usr/bin/env bash
# Build a symlink tree of calexp images for daomop-sns (header-based stacking).
#
# Usage:
#   daomop-link-calexps.sh LISTFILE OUTDIR [--basedir DIR] [--dry-run]
#
# LISTFILE lines are relative paths (one per line), e.g.:
#   CFHT_LSST/rerun/processCcdOutputs/calexp/23AP34/JF1/2023-02-24/gri/calexp-2847528-18.fits
# Field may contain spaces. Groups by (field, night, ccd); REFNUM = middle of sorted
# unique expnums (n//2).
# Links: OUTDIR/field/night/ccdXX/calexp-{EXPNUM}-{REFNUM}-{XX}.fits
#   ->   BASEDIR/<path_from_list>
set -euo pipefail

usage() {
    cat <<EOF
Usage: $(basename "$0") LISTFILE OUTDIR [--basedir DIR] [--dry-run]

  LISTFILE    List of relative calexp paths (one per line)
  OUTDIR      Root directory for symlink tree
  --basedir   Prefix for list paths (default: .)
  --dry-run   Print ln -s actions without creating links
EOF
    exit 1
}

BASEDIR="."
DRY_RUN=0
LISTFILE=""
OUTDIR=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --basedir)
            [[ $# -ge 2 ]] || usage
            BASEDIR="$2"
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
[[ -d "${BASEDIR}" ]] || { echo "BASEDIR not found: ${BASEDIR}" >&2; exit 1; }

# Resolve BASEDIR to absolute for symlink targets.
if command -v realpath >/dev/null 2>&1; then
    BASEDIR_ABS="$(realpath "${BASEDIR}")"
else
    BASEDIR_ABS="$(cd "${BASEDIR}" && pwd)"
fi

WORKDIR="$(mktemp -d "${TMPDIR:-/tmp}/daomop-link-calexps.XXXXXX")"
cleanup() { rm -rf "${WORKDIR}"; }
trap cleanup EXIT

# Parse list -> tab-separated: field, night, ccd, expnum, relpath
# Split on '/' only so field names with spaces stay intact.
awk '
{
    gsub(/\r/, "")
    sub(/^[[:space:]]+/, "")
    sub(/[[:space:]]+$/, "")
    if ($0 == "" || $0 ~ /^#/) next

    path = $0
    n = split(path, p, "/")
    if (n < 5) {
        print "Warning: skipping short path: " path > "/dev/stderr"
        next
    }
    base = p[n]
    if (base !~ /^calexp-[0-9]+-[0-9]{2}\.fits$/) {
        print "Warning: skipping unexpected basename: " base > "/dev/stderr"
        next
    }
    # calexp-EXPNUM-CCD.fits
    split(base, b, /[-.]/)
    # b[1]=calexp, b[2]=expnum, b[3]=ccd, b[4]=fits
    expnum = b[2]
    ccd = b[3]
    night = p[n - 2]
    field = p[n - 3]
    if (night !~ /^[0-9]{4}-[0-9]{2}-[0-9]{2}$/) {
        print "Warning: skipping path with unexpected night dir: " path > "/dev/stderr"
        next
    }
    printf "%s\t%s\t%s\t%s\t%s\n", field, night, ccd, expnum, path
}
' "${LISTFILE}" > "${WORKDIR}/raw.txt"

[[ -s "${WORKDIR}/raw.txt" ]] || {
    echo "No calexp entries found in ${LISTFILE}" >&2
    exit 1
}

# Unique by field,night,ccd,expnum (keep first path); then sort for grouping.
# Sort keys: field, night, ccd, numeric expnum
sort -t $'\t' -u -k1,1 -k2,2 -k3,3 -k4,4n "${WORKDIR}/raw.txt" > "${WORKDIR}/pairs.txt"

# Emit one record per (field, night, ccd):
#   field<TAB>night<TAB>ccd<TAB>REFNUM<TAB>expnum:path,expnum:path,...
awk -F '\t' '
{
    key = $1 "\t" $2 "\t" $3
    en = $4
    path = $5
    if (key != prev && prev != "") {
        n = length(exps)
        ref = exps[int(n / 2)]
        printf "%s\t%s", prev, ref
        for (i = 0; i < n; i++) {
            printf "\t%s:%s", exps[i], paths[exps[i]]
        }
        printf "\n"
        delete exps
        delete paths
    }
    if (!(en in paths)) {
        exps[length(exps)] = en
        paths[en] = path
    }
    prev = key
}
END {
    if (prev != "") {
        n = length(exps)
        ref = exps[int(n / 2)]
        printf "%s\t%s", prev, ref
        for (i = 0; i < n; i++) {
            printf "\t%s:%s", exps[i], paths[exps[i]]
        }
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

while IFS=$'\t' read -r field night ccd REFNUM rest; do
    # rest is tab-separated expnum:path entries; reconstruct from remaining fields
    # After REFNUM, remaining columns are expnum:path (path may contain no tabs)
    n_groups=$((n_groups + 1))

    # Count entries for the log line
    entry_count=0
    IFS=$'\t' read -r -a entries <<< "${rest}"
    entry_count=${#entries[@]}

    echo "Group ${field}/${night}/ccd${ccd}: ${entry_count} expnum(s), REFNUM=${REFNUM}"

    dest_dir="${OUTDIR}/${field}/${night}/ccd${ccd}"
    if [[ "${DRY_RUN}" -eq 0 ]]; then
        mkdir -p "${dest_dir}"
    else
        echo "mkdir -p ${dest_dir}"
    fi

    for entry in "${entries[@]}"; do
        [[ -z "${entry}" ]] && continue
        expnum="${entry%%:*}"
        relpath="${entry#*:}"
        src="${BASEDIR_ABS}/${relpath}"
        link_name="calexp-${expnum}-${REFNUM}-${ccd}.fits"
        dest="${dest_dir}/${link_name}"

        if [[ ! -e "${src}" ]]; then
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
done < "${WORKDIR}/groups.txt"

echo "----"
echo "Groups: ${n_groups}"
echo "Links created: ${n_links}"
echo "Missing sources: ${n_missing}"
if [[ "${DRY_RUN}" -eq 1 ]]; then
    echo "(dry-run: no files were created)"
fi
