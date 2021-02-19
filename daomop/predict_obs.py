import argparse
import os

import numpy
from astropy.table import Table
from astropy.time import Time
from mp_ephem import BKOrbit

from . import util


def run():
    parser = argparse.ArgumentParser()
    parser.add_argument('mpcfile', type=str)
    parser.add_argument('trackfile', type=str)
    parser.add_argument('--skycat-file', type=str, default='skycat.csv')
    args = parser.parse_args()
    main(**vars(args))


def main(mpcfile, trackfile, skycat_file):
    """

    :param mpcfile: Name of file containing MPC formated astrometric measurements.
    :param skycat_file: Name of file with table of field locations.
    :return:
    """
    skycat = Table.read(skycat_file, format='csv')
    required_names = ['mjdobs', 'ramin', 'ramax', 'decmin', 'decmax', 'pointing', 'ccd']
    for name in required_names:
        if name not in skycat.colnames:
            raise ValueError(f"{skycat_file} missing one or more columns: {required_names}")

    colnames = ["pointing",
                "ccd",
                "index",
                "x",
                "y",
                "rate",
                "angle",
                "visit",
                "stack",
                "ra",
                "dec",
                "nimg"]
    track_rows = []

    orbit = BKOrbit(None, ast_filename=mpcfile)
    name = orbit.observations[0].provisional_name
    try:
        pointing, ccd, index = util.from_provisional_name(name)
    except:
        index = name
    for row in skycat:
        orbit.predict(Time(row['mjdobs'], format='mjd'))
        coord1 = orbit.coordinate
        if not (row['ramin'] < coord1.ra.degree < row['ramax'] and
                row['decmin'] < coord1.dec.degree < row['decmax']):
            continue
        orbit.predict(Time(row['mjdobs']+1/24.0, format='mjd'))
        coord2 = orbit.coordinate
        ra_rate = (coord2.ra - coord1.ra).to('arcsec')
        dec_rate = (coord2.dec - coord1.dec).to('arcsec')
        angle = numpy.degrees(numpy.arctan2(dec_rate, ra_rate)).value
        if angle < 0:
            angle = angle + 180 % 360
        rate = (ra_rate**2 + dec_rate**2)**0.5
        track_rows.append([f"0{row['pointing']}",
                           row['ccd'],
                           index,
                           0,
                           0,
                           rate.value,
                           angle,
                           0,
                           0,
                           orbit.ra.degree,
                           orbit.dec.degree,
                           3])
    track_rows = numpy.array(track_rows)
    Table(track_rows, names=colnames).write(trackfile, format='ascii.commented_header')


if __name__ == '__main__':
    run()
