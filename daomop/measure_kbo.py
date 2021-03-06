"""
Create an Observation using ds9 displaying an image of a KBO source.
"""
import argparse
import logging
import os

import numpy
import pyds9
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from astropy.wcs import WCS
from mp_ephem import BKOrbit, EphemerisReader
from mp_ephem.ephem import Observation
from vos import Client

from . import settings
from . import util

config = settings.AppConfig()

field_ids = {3068: 'P68', 3071: 'P71', 3072: 'P72'}


def start_ds9(name):
    ds9 = pyds9.DS9(target=name)
    levels = ['INIT', 'PREF']
    for level in levels:
        setting = config.read(f"DS9.{level}")
        for key in list(setting.keys()):
            ds9.set(f"{key.replace('_', ' ')} {setting[key]}")
    ds9.set("frame delete all")
    return ds9


def get_ds9(name):
    return pyds9.DS9(target=name, start=False)


def main(orbit=None, **kwargs):
    """

    :param kwargs:
    :type orbit: BKOrbit
    :return:
    """
    from daomop import daophot
    pointing = kwargs['pointing']
    index = kwargs['index']
    ccd = kwargs['ccd']
    rate = kwargs['rate']
    angle = kwargs['angle']
    ra = kwargs['ra']
    dec = kwargs['dec']
    p_name = kwargs['p_name']
    discovery = kwargs['discovery']
    nimg = kwargs['nimg']
    # orbit = kwargs.get('orbit', None)
    # isinstance(BKOrbit, orbit)
    client = Client()

    int_rate = int(rate * 10)
    int_angle = int((angle % 360) * 10)
    images = []
    ds9 = get_ds9('validate')
    # Load the 3 images associated with this point/ccd/rate/angle set.
    for idx in range(nimg):
        expnum = f'{int(pointing)}{int_rate:02d}{int_angle:04d}{idx}'
        image = f'{expnum}p{ccd:02d}.fits'
        url = f'vos:NewHorizons/{index}/dbimages/{expnum}/ccd{ccd:02d}/{image}'
        try:
            if not os.access(image, os.R_OK):
                # get from VOSpace is not already on disk
                client.copy(url, image)
        except Exception as ex:
            logging.error(str(ex))
            # Return empty set on VOSpace copy error.
            return {}
        images.append(image)

    wcs_dict = {}
    ds9.set('frame delete all')
    ds9.set('zscale')
    for image in images:
        ds9.set('frame new')
        ds9.set(f"file {image}")
        with fits.open(image) as hdulist:
            header = hdulist[1].header
            obsdate = Time(hdulist[0].header['DATE-AVG'], scale='tai').utc
            wcs_dict[image] = WCS(header)
            if orbit is not None:
                orbit.predict(obsdate)
                ra = orbit.coordinate.ra.degree
                dec = orbit.coordinate.dec.degree
                uncertainty_ellipse = (orbit.dra.to('arcsec').value,
                                       orbit.ddec.to('arcsec').value,
                                       orbit.pa.to('degree').value + 90)
            else:
                uncertainty_ellipse = 3, 3, 0
            ds9.set('regions', f'icrs; ellipse({ra},{dec},'
                               f'{uncertainty_ellipse[0]}",'
                               f'{uncertainty_ellipse[1]}",'
                               f'{uncertainty_ellipse[2]})')
            ds9.set(f'pan to {ra} {dec} wcs icrs')
    ds9.set('frame match wcs')
    ds9.set('frame first')
    obs = {}
    # Build a map of allowed key strokes
    allowed_keys = {'x': ('', 'centroid at this location'),
                    'q': ('', 'Quit this image set'),
                    'p': ('', 'Previous frame'),
                    'n': ('', 'Next frame'),
                    'r': ('', 'Create a NULL observation.')}
    for key in [x.split() for x in config.read("MPC.NOTE1OPTIONS")]:
        allowed_keys[key[0].lower()] = key

    while True:
        try:
            result = ds9.get('iexam key coordinate image')
            key, x, y = result.split()
            x = float(x)
            y = float(y)
            logging.debug(f"DS9 Returned: {result} -> {key} {x} {y}")
        except Exception as ex:
            logging.debug(f"DS9 get exception: {ex}")
            continue

        if key == 'n':
            ds9.set('frame next')
            continue
        if key == 'p':
            ds9.set('frame previous')
            continue

        if key not in allowed_keys:
            logging.info(f"Allowed keys: ")
            for key in allowed_keys:
                print(f"{key} -> {allowed_keys[key][1]}")
            continue

        if key == 'q':
            break
        note1 = allowed_keys[key][0]
        frame_no = int(ds9.get('frame')) - 1
        image = images[frame_no]
        ds9.set('regions', f'image; circle {x} {y} 20')
        centroid = not note1 == 'H'
        phot = daophot.phot_mag(image,
                                [x, ], [y, ],
                                aperture=5,
                                sky_inner_radius=15,
                                sky_annulus_width=10,
                                apcor=0.3,
                                zmag=27.3,
                                maxcount=1000,
                                extno=1,
                                exptime=90.0,
                                centroid=centroid)

        phot_failure = (phot['PIER'][0] != 0 or
                        phot.mask[0]['MAG'] or
                        phot.mask[0]['MERR'])
        sky_failure = phot['SIER'][0] != 0
        cen_failure = phot['CIER'][0] != 0

        if phot_failure or sky_failure or cen_failure:
            logging.warning(f"iraf.daophot.phot error:\n {phot}")
            cen_x = x
            cen_y = y
            obs_mag = None
            obs_mag_err = None
            note1 = "H"
        else:
            cen_x = phot['XCENTER'][0]
            cen_y = phot['YCENTER'][0]
            obs_mag = phot['MAG'][0]
            obs_mag_err = phot['MERR'][0]

        obsdate = Time(Time(fits.open(image)[0].header['DATE-AVG'], scale='tai').mjd,
                       format='mjd',
                       precision=6).mpc
        try:
            ra, dec = wcs_dict[image].all_pix2world(cen_x, cen_y, 1)
        except Exception as ex:
            logging.warning(f"Failure converting {cen_x, cen_y} to RA/DEC for {image}")
            logging.warning(ex)
            logging.warning(f"Got: {ra},{dec}")

        record_key = obsdate[0:13]
        obs[record_key] = (Observation(
            null_observation=key == 'r',
            provisional_name=p_name,
            note1=note1,
            note2='C',
            date=obsdate,
            ra=ra,
            dec=dec,
            mag=obs_mag,
            mag_err=obs_mag_err,
            band='r',
            observatory_code='568',
            discovery=discovery,
            comment=None,
            xpos=x,
            ypos=y,
            frame=image,
            astrometric_level=0))
        discovery = False
    return obs


def _main(**kwargs):
    start_ds9('validate')
    if 'provisional_name' not in kwargs:
    #    kwargs['provisional_name'] = util.get_provisional_name(**kwargs)
         kwargs['provisional_name'] = kwargs['index']
    ast_filename = f"{kwargs['provisional_name']}.mpc"
    logging.info(f"Attempting measures of {kwargs['provisional_name']}, will write to {ast_filename}")
    obs = {}
    kwargs['discovery'] = True
    if os.access(ast_filename, os.R_OK):
        if kwargs['skip']:
            logging.warning(f"{ast_filename} exists, skipping")
            return
        logging.warning(f"{ast_filename} exists, appending new measures.")
        for ob in list(EphemerisReader().read(ast_filename)):
            record_index = ob.date.mpc[0:13]
            obs[record_index] = ob
        try:
            kwargs['orbit'] = BKOrbit(None, ast_filename=ast_filename)
            kwargs['discovery'] = False
        except Exception as er:
            logging.warning(f"{ast_filename} -> {er}")

    new_obs = main(**kwargs)
    logging.info(f"{new_obs}")
    for record_index in new_obs:
        if record_index in obs:
            # preserve the discovery key
            if obs[record_index].discovery:
                new_obs[record_index].discovery = True
        for old_record_index in obs:
            # remove measures in obs that have the same frame as the new obs.
            if new_obs[record_index].comment.frame == obs[old_record_index].comment.frame:
                del obs[old_record_index]
        obs[record_index] = new_obs[record_index]

    obs_list = [obs[record_index] for record_index in obs]
    with open(ast_filename, 'w') as mpc_obj:
        for ob in obs_list:
            mpc_obj.write(ob.to_string() + "\n")
    try:
        orbit = BKOrbit(None, ast_filename=ast_filename)
        logging.info(orbit.summarize())
    except Exception as ex:
        logging.error(f"{ex}")


def get_valid_obs_count(observations):
    nobs = 0
    for obs in observations:
        if not obs.null_observation:
            nobs += 1
    return nobs


def main_args(args):
    kwargs = vars(args)
    kwargs['p_name'] = kwargs['index']
    # kwargs['p_name'] = util.get_provisional_name(**kwargs)
    _main(**kwargs)


def main_list(args):
    t = Table.read(args.detections, format='ascii')
    logging.debug(f"read table of len {len(t)} with columns {t.colnames}")

    targets = numpy.unique(t['index'])
    if len(targets) == 1:
        t['provisional_name'] = os.path.splitext(os.path.basename(args.detections))[0]
        logging.info(f"Using provisional name {t['provisional_name'][0]}")

    for target in targets:
        rows = t[t['index'] == target]

        if 'provisional_name' not in rows.colnames:
            logging.debug(f"Computing provisional name for {rows[0]}")
            kwarg = dict([(name, rows[0][name]) for name in rows[0].colnames])
            p_name = util.get_provisional_name(**kwarg)
            logging.debug(f"Got: {p_name}")
        else:
            p_name = rows[0]['provisional_name']

        logging.info(f"{p_name} available on {len(rows)} nights")

        for row in rows:
            logging.debug(f"{row}")
            kwarg = dict([(name, row[name]) for name in row.colnames])
            kwarg['p_name'] = p_name
            kwarg['skip'] = args.skip
            logging.debug(kwarg)
            _main(**kwarg)


def run():
    main_parser = argparse.ArgumentParser()
    subparsers = main_parser.add_subparsers(help='sub-command help')
    parser_file = subparsers.add_parser('file', help='a help')
    parser_file.add_argument('detections', help='Name of file containing list of detections, must have columns'
                                                'for pointing, ccd, index, rate, angle, ra, dec')
    parser_file.set_defaults(func=main_list)
    parser_args = subparsers.add_parser('args', help='a help')
    parser_args.add_argument('pointing', type=int)
    parser_args.add_argument('ccd', type=int)
    parser_args.add_argument('index', help="index of the detection" )
    parser_args.add_argument('x', type=float)
    parser_args.add_argument('y', type=float)
    parser_args.add_argument('rate', type=float)
    parser_args.add_argument('angle', type=float)
    parser_args.add_argument('visit', type=int)
    parser_args.add_argument('stack', type=str)
    parser_args.add_argument('ra', type=float, help="RA of initial object location (deg)")
    parser_args.add_argument('dec', type=float, help="DEC of initial object location (deg)")
    parser_args.add_argument('nimg', type=int, help="Number of images in stack set")
    parser_args.set_defaults(func=main_args)
    main_parser.add_argument('--log-level', choices=['DEBUG', 'INFO', 'ERROR'], default='INFO')
    main_parser.add_argument('--skip', action="store_true")

    args = main_parser.parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level))
    args.func(args)


if __name__ == '__main__':
    run()
