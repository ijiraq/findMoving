import argparse
import logging
import string

import numpy as np
import re
from pathlib import Path
import os

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


def get_provisional_name(pointing, ccd, index, **kwargs):
    """
    Get a provisional name based on the pointing, ccd and index values

    value is Letter for the pointing '3000' ==> P '3100' ==> Q etc.
             and HEX encoded value of {ccd:03d}{index:04d}

    e.g.  3091 67 1 => P97a3931

    :param pointing:  LSST pointing number
    :param ccd: CCD on which source is detected
    :param index: INDEX of detection on that CCD
    :param kwargs: Not used
    :return:
    """
    logging.debug(f"{pointing} {ccd} {index} {kwargs}")
    try:
        index = int(index)
    except:
        index = 1
    L = (string.digits + string.ascii_uppercase)[pointing//100-5]
    value = int(f"{ccd:03d}{index:04d}")
    return f"{L}{str(pointing)[-2:]:2s}{hex(value)[2:]:4s}"


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
        _value = str(int(f'0x{p_name[-5:]}', base=16))
        ccd = _value[0:2]
        index = _value[2:2+4]

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
