import argparse
import logging
import string

import numpy as np
import re
from pathlib import Path
import os
from astropy import time
from astropy.time import Time

base_parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                      fromfile_prefix_chars='@',
                                      add_help=False)
base_parser.add_argument('basedir', help="Root directory of LSST pipelined data")
base_parser.add_argument('--pointing', help="Which pointing to process, eg. 03071", default=None)
base_parser.add_argument('--field', help='Which FIELD to process', choices=['NHF1', 'NHF2'], default=None)
base_parser.add_argument('--patch', help="sky_inner_radius patch to process (e.g. 0,0)", nargs=1)
base_parser.add_argument('--exptype', help="PREFIX to select type of exposure",
                         type=str, default='DIFFEXP')
base_parser.add_argument('--rerun', help="rerun directory containing the warped difference images.", nargs=1)
base_parser.add_argument('--filter', help="Filter to stack", default="HSC-R2")
base_parser.add_argument('--visit', type=int, help='visit to process.', default=None)
base_parser.add_argument('--ccd', help="Which CCD to stack?", type=int, default=None)
base_parser.add_argument('--log-level', help="What level to log at?", default="ERROR",
                         choices=['INFO', 'ERROR', 'DEBUG'])


def get_image_list(dirname, exptype='CORR', visit=None, ccd=None, filters=[]):
    _filelist = []
    exptype = exptype is not None and f"{exptype}" or "*"
    visit = visit is not None and f"{visit:07d}" or "???????"
    ccd = ccd is not None and f"{ccd:03d}" or "???"
    pattern = f"{exptype}-{visit}-{ccd}.f*"
    logging.info(f"Searching for data in {dirname} -> {Path(dirname).resolve()} using pattern: {pattern}")
    for path in Path(dirname).rglob(f'{pattern}'):
        if not re.search(exptype+'-[0-9]{7}-[0-9]{3}.fits', path.name):
            continue
        full_path = str(path.resolve())
        filtered_in = True
        for pattern in filters:
            logging.debug(f"Filtering {full_path} using pattern: {pattern}")
            if not re.search(str(pattern), full_path):
                filtered_in = False
                break
        if filtered_in:
            _filelist.append(full_path)
    _filelist = np.unique(_filelist)
    logging.info(f"Returning list of {len(_filelist)} files.")
    return _filelist


def get_provisional_name(pointing, chip, index, **kwargs):
    """
    Get a provisional name based on the pointing, ccd and index values

    value is Letter for the pointing '3000' ==> P '3100' ==> Q etc.
             and HEX encoded value of {ccd:03d}{index:04d}

    e.g.  3091 67 1 => P97a3931

    :param pointing:  LSST pointing number
    :param chip: CCD on which source is detected
    :param index: INDEX of detection on that CCD
    :param kwargs: Not used
    :return:
    """
    logging.debug(f"{pointing} {chip} {index} {kwargs}")
    try:
        index = int(index)
    except:
        index = 1
    try:
        ccd = int(chip)
    except:
        ccd = 0
    L = (string.digits + string.ascii_uppercase)[pointing//100-5]
    value = (hex(int(f"{ccd:03d}{index:02d}"))[2:]).rjust(4, '0')  # 10301 --> 283d
    return f"{L}{str(pointing)[-2:]:2s}{value:4s}"


def from_provisional_name(p_name):
    """
    Get the pointing, ccd and index from the provisional name

    P97a3931 => 3097 67 1   (a3931) -> 0xa3931 -> turned in an integer -> first 3 digits are CCD next 4 are INDEX
    P => 3000
    Q => 3100
    etc.
    :param p_name:
    :return: pointing, ccd, index
    """
    if p_name[3] == 'x':
        index = int(f'0{p_name[3:]}', base=16)
        ccd = None
    else:
        _value = str(int(f'0x{p_name[3:]}', base=16)).rjust(5, '0')
        ccd = int(_value[0:3])
        index = int(_value[3:])
        if ccd > 103:
            ccd = int(_value[0:2])
            index = int(_value[2:])

    pointing = ((string.digits + string.ascii_uppercase).find(p_name[0])+5) * 100 + int(p_name[1:3])

    return pointing, ccd, index


def parse_rerun(basedir, rerun):
    reruns = rerun[0].split(":")
    if len(reruns) > 2:
        raise ValueError("Don't know what to do with more then 2 rerun directories.")

    if len(reruns) == 1:
        input_rerun = output_rerun = reruns[0]
    else:
        input_rerun = reruns[0]
        output_rerun = reruns[1]
    input_rerun = os.path.join(basedir, 'rerun', input_rerun)
    output_rerun = os.path.join(basedir, 'rerun', output_rerun)

    return input_rerun, output_rerun


def mid_exposure_mjd(hdu):
    header = hdu.header
    if 'DATE-AVG' in header:
        return time.Time(header['DATE-AVG'], scale='tai')
    mjd_start = time.Time(header['MJD-OBS'], format='mjd', scale='utc')
    mjd_end = time.Time(header['MJDEND'], format='mjd', scale='utc')
    return (mjd_start + (mjd_end - mjd_start) / 2.0).tai


def mid_exposure_mpc(hdu, precision=5):
    return Time(mid_exposure_mjd(hdu).utc, format='mjd', precision=precision).mpc


def exposure_id_from_name(name):
    """
    Parse exposure id from an image basename.

    Supports MegaCam/CFHT diff and calexp-style names, e.g.:
      DIFFEXP-2779494-00.fits
      DIFFEXP-2779494-2779500-00.fits
      calexp-02779494-012.fits
    """
    base = os.path.basename(str(name))
    match = re.search(
        r'(?:DIFFEXP|calexp|CORR|WARPED)[-_]?(\d{6,8})',
        base,
        flags=re.IGNORECASE,
    )
    if match:
        return int(match.group(1))
    for part in Path(base).stem.split('-'):
        if part.isdigit() and len(part) >= 6:
            return int(part)
    return None


def exposure_id_from_filename_keyword(filename):
    """
    Convert a CFHT/MegaCam FILENAME keyword (e.g. ``2779465o``) to an exposure id.
    """
    if filename is None or filename == '':
        return None
    value = str(filename).strip()
    if value.lower().endswith('o'):
        value = value[:-1]
    try:
        return int(value)
    except ValueError:
        return value or None


def frame_label(hdu, fallback_name=None):
    """
    MegaCam exposure id for logging / headers.

    Prefers header ``EXPNUM``; if missing, uses ``FILENAME`` with a trailing
    ``o`` removed (e.g. ``2779465o`` → ``2779465``). Falls back to parsing
    ``IMAGE`` / ``fallback_name`` only if those keywords are absent.
    """
    try:
        header = hdu[0].header
    except (TypeError, KeyError, IndexError, AttributeError):
        header = hdu.header

    if 'EXPNUM' in header and header['EXPNUM'] not in (None, ''):
        try:
            return int(header['EXPNUM'])
        except (TypeError, ValueError):
            return header['EXPNUM']

    if 'FILENAME' in header and header['FILENAME'] not in (None, ''):
        exp_id = exposure_id_from_filename_keyword(header['FILENAME'])
        if exp_id is not None:
            return exp_id

    name = header.get('IMAGE') or fallback_name or ''
    exp_id = exposure_id_from_name(name)
    if exp_id is not None:
        return exp_id
    return os.path.basename(str(name)) if name else 'unknown'


_NIGHT_DIR_PATTERN = re.compile(r'^\d{4}-\d{2}-\d{2}$')


def night_from_path(path):
    """
    Parse UT night (``YYYY-MM-DD``) from a directory path.

    Matches the field/night/ccdXX symlink layout used by daomop-link-diffs.sh.
    """
    if not path:
        return None
    for part in Path(str(path)).parts:
        if _NIGHT_DIR_PATTERN.match(part):
            return part
    return None


def night_from_hdu(hdu):
    """UT night of mid-exposure from a FITS primary HDU."""
    return mid_exposure_mjd(hdu).utc.strftime('%Y-%m-%d')


def stack_night(hdu, fallback_name=None):
    """
    UT night for a stack, from WARPD path or reference exposure time.
    """
    try:
        header = hdu[0].header
    except (TypeError, KeyError, IndexError, AttributeError):
        header = hdu.header

    night = night_from_path(header.get('WARPD'))
    if night is not None:
        return night
    if fallback_name:
        night = night_from_path(os.path.dirname(fallback_name))
        if night is not None:
            return night
    try:
        return night_from_hdu(hdu[0])
    except (TypeError, KeyError, IndexError, AttributeError):
        return night_from_hdu(hdu)


def field_from_path(path):
    """
    Parse field/pointing name from a directory path.

    Matches the field/night/ccdXX symlink layout (parent of the night folder).
    """
    if not path:
        return None
    parts = list(Path(str(path)).parts)
    for index, part in enumerate(parts):
        if _NIGHT_DIR_PATTERN.match(part) and index > 0:
            return parts[index - 1]
    return None


def sanitize_filename_component(name):
    """Make a string safe for use as one component of a stack output filename."""
    if name is None or name == '':
        return None
    value = str(name).strip()
    value = re.sub(r'[^\w\-.]+', '_', value)
    value = re.sub(r'_+', '_', value).strip('_')
    return value or None


def source_name_from_orbit(orbit):
    """
    Provisional designation from astrometry loaded into a BKOrbit.

    Uses ``orbit.observations[0].provisional_name``, falling back to ``orbit.name``
    if set during the orbit fit.
    """
    try:
        observations = orbit.observations
        if observations:
            name = getattr(observations[0], 'provisional_name', None)
            if name is not None and str(name).strip():
                return str(name).strip()
    except (AttributeError, IndexError, TypeError):
        pass
    name = getattr(orbit, 'name', None)
    if name is not None and str(name).strip():
        return str(name).strip()
    return None


def stack_output_filename(pointing, night, ccd, n_sub_stacks, index, stack_mode, rate, angle,
                          source_name=None):
    """
    Build a stack output filename aligned with field/night/ccd directory layout.

    Format:
      STACK-{source}-{pointing}-{night}-{ccd}-{n_sub_stacks}-{index}-{stack_mode}{rate}{angle}.fits
    ``source`` is omitted when ``source_name`` is not given.
    """
    field = sanitize_filename_component(pointing) or 'UNKNOWN'
    night_label = night or 'UNKNOWN'
    parts = ['STACK']
    source = sanitize_filename_component(source_name)
    if source:
        parts.append(source)
    parts.extend([
        field,
        night_label,
        f'{ccd:02d}',
        f'{n_sub_stacks:02d}',
        f'{index:02d}',
        f'{stack_mode}{rate:+06.2f}-{angle:+06.2f}',
    ])
    return '-'.join(parts) + '.fits'
