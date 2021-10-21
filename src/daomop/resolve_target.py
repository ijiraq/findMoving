"""
Look up the location of source in a stack given the x/y/rate/angle of the source and stack.

Summary:  using lines like this

chip x y mag ra_rate dec_rate
0.0 1362.14 3656.88 25.37 -1.74 -0.26

produce lines like this

 pointing  chip  index         x         y      rate     angle                                       stack          ra         dec  nstk
    03147   000    001   1362.14   3656.88      1.76     -8.50  STACK-223348-000-FU-+05.00--10.00.fits.fz    288.07598   -20.78413     3

The output can be used as input into daomop-target-sns.sh

"""
import argparse
import logging
import math
import os
import re
from io import BytesIO

import numpy
import vos
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

OUTPUT_COLUMNS = ["pointing", "chip", "index",   "x",   "y", "rate", "angle",  "mag", "stack", "ra",  "dec", "nstk"]
OUTPUT_COLUMN_FMT = {'pointing': '05d',
                     'chip': '03d',
                     'index': '03d',
                     'x': '8.2f',
                     'y': '8.2f',
                     'rate': '8.2f',
                     'angle': '8.2f',
                     'mag': '8.2f',
                     'stack': '>50s',
                     'ra': '10.5f',
                     'dec': '10.5f',
                     'nstk': '3d'}

RESOLVED_COLUMN_NAMES = ['stack', 'ra', 'dec', 'nstk']
RESOLVED_COLUMN_FILL = {'stack': ' '*50, 'ra': 0., 'dec': 0., 'nstk': 0}

CUTOUT = "[1][1:1,1:1]"


STACK_PATTERN = re.compile(
    r"STACK-(?P<visit>\d{6,7})-(?P<chip>\d{3})(_masked)?-00-(?P<rate>[+-]\d{2}.\d{2})-(?P<angle>[+-]\d{2}.\d{2}).fits.fz")


def resolve(pointing, chip, rate, angle, x, y, vos_basedir):
    """

    :param pointing: pointing that source is found on
    :param chip: chip that source is found on
    :param rate: rate of motion that maximized source signal
    :param angle: angle of motion that maximized source signal
    :param x: x pixel value of source location in pointing/chip/rate/angle stack
    :param y: y pixel value
    :param vos_basedir: base location of VOSpace stacks, e.g. 'vos:NewHorizons/S20A-OT04/STACKS_V4'
    :return: stack, ra, dec, mjd
    """

    # Resolve retrieves the headers from the VOSpace repo with the stacks.
    client = vos.Client()
    # Standard pattern used for DAOMOP Shift-and-Stack

    # Get a list of stack files and match to the one whose rate is closest to what we need.
    stack_dir_uri = f"{vos_basedir}/{pointing:05d}/{chip:03d}"
    stacks = client.listdir(stack_dir_uri)
    stack = None
    min_offset = None
    min_idx = None
    for idx, stack in enumerate(stacks):
        stack_angle = stack_rate = 0
        try:
            match = STACK_PATTERN.search(stack)
            stack_rate = float(match.group('rate'))
            stack_angle = float(match.group('angle'))
        except Exception as ex:
            # logging.debug(f"Failed to match {stack}\n")
            # logging.debug(str(ex))
            continue
        # use cos/sin rule to compute vector difference between requested rate and stack rate.
        da = math.radians(math.fabs(stack_angle-angle))
        dr = rate**2 + stack_rate**2 - 2*rate*stack_rate*math.cos(da)
        # select the minimum vector difference as the stack of interest.
        if min_offset is None or min_offset > dr:
            min_offset = dr
            min_idx = idx

    # Found the STACK in the VOSpace storage area that matches the pointing/chip/rate/angle of interest.
    # Now pull the 'cutout' with the WCS and compute the RA/DEC from the x/y location.
    stack_image = stacks[min_idx]
    logging.info(f"pointing:{pointing:05d}, chip:{chip:03d}, rate:{rate:5.1f}, angle:{angle:5.1f} "
                 f"matched stack:{stack_image} in {vos_basedir}")
    uri = f"{stack_dir_uri}/{stack_image}"

    # noinspection PyBroadException
    # get the header of the matching stack and compute RA/DEC of the source using that header
    try:
        fobj = BytesIO(client.open(uri, cutout=CUTOUT, view='cutout').read())
        with fits.open(fobj) as hdu:
            ra, dec = WCS(hdu[0].header).all_pix2world(x, y, 1)
    except Exception as ex:
        logging.error("Failed while measuring {pointing} {chip} {x} {y} on {stack}")
        logging.error("{ex}")
        return None

    return {'stack': stack_image, 'ra': ra, 'dec': dec}


def main():
    """
    The executable program, parses the args and calls resovle
    :return:
    """

    parser = argparse.ArgumentParser(description="""Compute the RA/DEC of candidate KBOs detected by WTF ML process.
    The input file is modified in place.""",
                                     epilog="detection file is a an ascii file with at least the following columns:"
                                            "'chip', 'id', 'x', 'y', 'rate', 'angle'")
    parser.add_argument('detection_filename', help="Name of detection list file from ML code.")
    parser.add_argument('--pointing', help="Sometimes the pointing isn't given in the input file.")
    parser.add_argument('--vospace-base-dir', default="vos:NewHorizons/S20A-OT04")
    parser.add_argument('--stack-version', default=4, help="Which version of the stacks are we looking at?")
    parser.add_argument('--nstack', help="how many cutout group stacks will be requested?", default=3)
    parser.add_argument('--log-level', default='INFO')

    args = parser.parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level))
    input_filename = args.detection_filename
    output_filename = f"{os.path.splitext(os.path.basename(input_filename))[0]}_out.txt"
    vospace_basedir = f"{args.vospace_base_dir}/STACKS_V{args.stack_version}"

    detections = Table.read(input_filename, format='ascii')
    if 'ra_rate' in detections.colnames and 'rate' not in detections.colnames:
        # convert ra_rate and dec_rate from planting process to expected values for stacking.
        detections['rate'] = (detections['ra_rate']**2 + detections['dec_rate']**2)**0.5
        detections['angle'] = numpy.degrees(numpy.arctan2(-1*detections['dec_rate'], numpy.fabs(detections['ra_rate'])))
    if 'index' not in detections.colnames:
        detections['index'] = -1

    # Some files record the chip as a floating point.
    detections['chip'] = numpy.int32(detections['chip'])

    if 'pointing' not in detections.colnames:
        detections['pointing'] = int(args.pointing)

    for colname in RESOLVED_COLUMN_NAMES:
        detections[colname] = RESOLVED_COLUMN_FILL[colname]
    idx = {}
    for row in detections:
        if 'NEW' in row and row['NEW'] == 'n':
            # this isn't a NEW source in this table.
            continue
        result = resolve(row['pointing'], row['chip'], row['rate'], row['angle'], row['x'], row['y'], vos_basedir=vospace_basedir)
        if row['index'] < 1:
            if row['chip'] not in idx:
                idx[row['chip']] = 0
            idx[row['chip']] += 1
            row['index'] = idx[row['chip']]

        if result is None:
            continue
        for colname in result:
            row[colname] = result[colname]

    detections['nstk'] = 3
    for colname in OUTPUT_COLUMN_FMT:
        detections[colname].info.format = OUTPUT_COLUMN_FMT[colname]
    detections[OUTPUT_COLUMNS].write(output_filename, format='ascii.fixed_width', delimiter=None, overwrite=True)
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())
