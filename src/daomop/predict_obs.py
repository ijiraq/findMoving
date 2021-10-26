"""
Given an ephemeris file and list of RA/DEC bounding boxes determine which boxes will contain the object on orbit determine from ephem.
"""
import argparse
import numpy
from astropy.table import Table
from astropy.time import Time
from mp_ephem import BKOrbit
import logging
from daomop import util


def run():
    """
    Execution entry point when running as script.
    :return:
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('mpc_filename', type=str, help="MPC formatted ephemeris file, will be passed to BKOrbit")
    parser.add_argument('track_filename', type=str, help="Nane of file to store lines that indicate where the source might be found.")
    parser.add_argument('--pointing-catalog-filename', type=str, default='skycat.csv',
                        help="File containing the bounding boxes for each chip of each pointing.")
    parser.add_argument('--log-level', default='ERROR', choices=['INFO', 'DEBUG', 'ERROR'])
    args = parser.parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level))
    main(**vars(args))


# noinspection PyBroadException
def main(mpc_filename, track_filename, pointing_catalog_filename, **kwargs):
    """

    :param mpc_filename: Name of file containing MPC formated astrometric measurements.
    :param track_filename: Name of file to write the tracks to.
    :param pointing_catalog_filename: Name of file with table of field locations.
    :return:
    """
    skycat = Table.read(pointing_catalog_filename, format='csv')
    required_names = ['mjdobs', 'ramin', 'ramax', 'decmin', 'decmax', 'pointing', 'ccd']
    for name in required_names:
        if name not in skycat.colnames:
            raise ValueError(f"{pointing_catalog_filename} missing one or more columns: {required_names}")

    colnames = ["pointing",
                "chip",
                "index",
                "x",
                "y",
                "rate",
                "angle",
                "mag",
                "stack",
                "ra",
                "dec",
                "nimg"]
    track_rows = []

    orbit = BKOrbit(None, ast_filename=mpc_filename)
    name = orbit.observations[0].provisional_name
    try:
        # use the logic of the provisional name to get the index
        pointing, ccd, index = util.from_provisional_name(name)
    except Exception as ex:
        logging.error(f"{ex} on {name}")
        index = 1
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
        angle = numpy.degrees(numpy.arctan2(dec_rate, numpy.fabs(ra_rate))).value
        rate = (ra_rate**2 + dec_rate**2)**0.5
        track_rows.append([row['pointing'],
                           row['ccd'],
                           index,
                           0.00,
                           0.00,
                           rate.value,
                           angle,
                           0.00,
                           "UNKNOWN",
                           orbit.ra.degree,
                           orbit.dec.degree,
                           3])
    track_rows = numpy.array(track_rows)
    Table(track_rows, names=colnames).write(track_filename, format='ascii.fixed_width', delimiter=None)


if __name__ == '__main__':
    run()
