"""
Query CADC SSOIS to find which MegaCam CHIP (ccdXX) directories contain a target.

Given an MPC/astrometry file and a search date window, calls the SSOIS arc search
(Bernstein fit) with extension resolution and reports unique ccd directories plus the
full TSV hit table.

SSOIS Ext is MegaCam FITS HDU index (1-based); ccd number is Ext - 1
(e.g. Ext=1 -> ccd00), matching cutouts like 2778831p00.fits.fz.
"""
from __future__ import annotations

import argparse
import logging
import sys
from io import StringIO
from typing import Iterable, Optional
from urllib.parse import urlencode

import requests
from astropy.io import ascii
from astropy.table import Table
from astropy.time import Time
from mp_ephem import EphemerisReader

SSOS_URL = "https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/cadcbin/ssos/ssosclf.pl"
DEFAULT_TELINST = "CFHT/MegaCam"
NEW_LINE = "\r\n"


def _parse_time(value: str) -> Time:
    """Parse a user epoch string as UTC Time (date or datetime)."""
    return Time(value, scale="utc")


def _epoch_for_ssois(value: str) -> str:
    """Convert a date/datetime string to SSOIS epoch form 'YYYY MM DD'."""
    t = _parse_time(value).replicate(format="iso")
    t.out_subfmt = "date"
    # iso date is YYYY-MM-DD; SSOIS wants spaces
    return str(t).replace("-", " ")


def load_observations(astfile: str) -> list:
    """Load MPC/ObsRecord lines from an astrometry file."""
    observations = EphemerisReader().read(astfile)
    if observations is None or len(observations) == 0:
        raise ValueError(f"No observations loaded from {astfile}")
    # Drop null observations if present
    cleaned = []
    for obs in observations:
        if getattr(obs, "null_observation", False):
            continue
        cleaned.append(obs)
    if not cleaned:
        raise ValueError(f"All observations in {astfile} are null")
    return cleaned


def query_ssois(
    observations: Iterable,
    epoch1: str,
    epoch2: str,
    telinst: str = DEFAULT_TELINST,
    search: str = "bern",
    extres: bool = True,
    xyres: bool = False,
    timeout: float = 180.0,
) -> Table:
    """
    Query SSOIS and return a Table of hits.

    Uses GET against ssosclf.pl with format=tsv (mp_ephem.ssos still points at a
    retired ssos.pl endpoint).
    """
    obs_blob = NEW_LINE.join(str(obs) for obs in observations)
    params = {
        "lang": "en",
        "obs": obs_blob,
        "search": search,
        "epoch1": _epoch_for_ssois(epoch1),
        "epoch2": _epoch_for_ssois(epoch2),
        "eunits": "none",
        "extres": "yes" if extres else "no",
        "xyres": "yes" if xyres else "no",
        "optical": "on",
        "ir": "on",
        "submm": "on",
        "telinst": telinst,
        "format": "tsv",
    }
    logging.info("SSOIS query %s epoch1=%s epoch2=%s telinst=%s",
                 SSOS_URL, params["epoch1"], params["epoch2"], telinst)
    logging.debug("SSOIS URL with params: %s?%s", SSOS_URL, urlencode(params))

    response = requests.get(SSOS_URL, params=params, timeout=timeout,
                            headers={"User-Agent": "daomop-ssos-chips"})
    if response.status_code != requests.codes.ok:
        raise IOError(f"SSOIS HTTP {response.status_code}: {response.text[:500]}")

    text = response.content.decode("utf-8", errors="replace")
    if "An error occured getting the ephemeris" in text:
        raise IOError(f"SSOIS ephemeris error:\n{text[:1000]}")

    # Skip any diagnostic preamble; table starts at Image\t header.
    lines = text.splitlines()
    start = 0
    for i, line in enumerate(lines):
        if line.startswith("Image\t") or line.startswith("Image "):
            start = i
            break
    else:
        raise IOError(f"SSOIS response missing Image header:\n{text[:1000]}")

    table_text = "\n".join(lines[start:])
    table = ascii.read(table_text.splitlines(), format="tab", guess=False, fast_reader=False)
    return table


def ext_to_ccd(ext) -> Optional[str]:
    """
    Map SSOIS Ext (FITS HDU index) to zero-padded ccd id.

    Returns None if Ext is missing / non-numeric (whole-mosaic hit).
    """
    if ext is None:
        return None
    s = str(ext).strip()
    if s in ("", "None", "nan", "--"):
        return None
    try:
        ext_i = int(float(s))
    except (TypeError, ValueError):
        return None
    if ext_i < 1:
        return None
    return f"{ext_i - 1:02d}"


def annotate_chips(table: Table) -> Table:
    """Add a 'ccd' column (zero-padded string) derived from Ext."""
    ccds = [ext_to_ccd(ext) for ext in table["Ext"]]
    out = table.copy()
    out["ccd"] = ccds
    return out


def filter_by_epoch(table: Table, epoch1: str, epoch2: str) -> Table:
    """
    Keep only rows whose MJD lies within [epoch1, epoch2].

    SSOIS accepts calendar dates only and often returns a wider night range
    than requested; apply the caller's full datetime window here.
    """
    if len(table) == 0 or "MJD" not in table.colnames:
        return table

    t1 = _parse_time(epoch1)
    t2 = _parse_time(epoch2)
    if t2 < t1:
        raise ValueError(f"epoch2 ({epoch2}) is before epoch1 ({epoch1})")

    mjd = table["MJD"].astype(float)
    mask = (mjd >= t1.mjd) & (mjd <= t2.mjd)
    n_in = int(mask.sum())
    logging.info(
        "Filtered SSOIS hits to MJD [%.5f, %.5f] (%s .. %s): %d of %d kept",
        t1.mjd, t2.mjd, t1.iso, t2.iso, n_in, len(table),
    )
    return table[mask]


def unique_ccd_dirs(table: Table) -> list[str]:
    """Sorted unique ccdXX directory names from an annotated table."""
    found = set()
    for ccd in table["ccd"]:
        if ccd is None or str(ccd) in ("None", "nan"):
            continue
        found.add(f"ccd{ccd}")
    return sorted(found)


def unique_night_ccd_dirs(table: Table) -> list[str]:
    """
    Sorted unique Image_target/yyyy-mm-dd/ccdXX paths from the hit table.

    Matches the field/night/ccdXX portion of the symlink tree layout.
    """
    if len(table) == 0:
        return []
    found = set()
    for row in table:
        ccd = row["ccd"]
        if ccd is None or str(ccd) in ("None", "nan"):
            continue
        night = Time(float(row["MJD"]), format="mjd", scale="utc").iso[:10]
        field = str(row["Image_target"]).strip()
        if field in ("", "None", "nan"):
            field = "UNKNOWN"
        found.add(f"{field}/{night}/ccd{ccd}")
    return sorted(found)


def main(
    astfile: str,
    epoch1: str,
    epoch2: str,
    telinst: str = DEFAULT_TELINST,
    output: Optional[str] = None,
    chips_only: bool = False,
    **_kwargs,
) -> int:
    observations = load_observations(astfile)
    provisional = getattr(observations[0], "provisional_name", None) or str(observations[0])[:12].strip()
    table = query_ssois(observations, epoch1, epoch2, telinst=telinst)
    table = filter_by_epoch(table, epoch1, epoch2)
    table = annotate_chips(table)
    night_chips = unique_night_ccd_dirs(table)

    logging.info("%s: %d SSOIS hit(s) in window, paths: %s",
                 provisional, len(table), ", ".join(night_chips) or "(none)")

    if chips_only:
        text = "\n".join(night_chips) + ("\n" if night_chips else "")
    else:
        # Prefer a compact summary then the table
        header = (
            f"# astfile={astfile} provisional={provisional}\n"
            f"# epoch1={epoch1} epoch2={epoch2} telinst={telinst}\n"
            f"# chips={' '.join(night_chips) if night_chips else '(none)'}\n"
            f"# nhits={len(table)} (SSOIS date search post-filtered by MJD window)\n"
        )
        buf = StringIO()
        table.write(buf, format="ascii.tab")
        text = header + buf.getvalue()

    if output:
        with open(output, "w", encoding="utf-8") as fh:
            fh.write(text)
        logging.info("Wrote %s", output)
    else:
        sys.stdout.write(text)

    return 0


def run():
    parser = argparse.ArgumentParser(
        description="Find MegaCam ccdXX directories for a target via CADC SSOIS.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("astfile", help="MPC / tnodb astrometry file for the target")
    parser.add_argument("--epoch1", required=True,
                        help="Search window start (UTC date or datetime, e.g. 2022-08-22T04:00:00)")
    parser.add_argument("--epoch2", required=True,
                        help="Search window end (UTC date or datetime, e.g. 2022-08-22T16:00:00)")
    parser.add_argument("--telinst", default=DEFAULT_TELINST,
                        help="SSOIS telescope/instrument filter")
    parser.add_argument("-o", "--output", default=None,
                        help="Write results to this file instead of stdout")
    parser.add_argument("--chips-only", action="store_true",
                        help="Print only unique Image_target/yyyy-mm-dd/ccdXX paths (one per line)")
    parser.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "ERROR"])
    args = parser.parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level), format="%(message)s")
    try:
        raise SystemExit(main(**vars(args)))
    except Exception as exc:
        logging.error("%s", exc)
        raise SystemExit(1) from exc


if __name__ == "__main__":
    run()
