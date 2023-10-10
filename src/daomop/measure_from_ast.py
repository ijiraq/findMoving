"""
Create an Observation using ds9 displaying an image of a KBO source.
"""
from mp_ephem import BKOrbit, EphemerisReader
from mp_ephem.ephem import Observation
import argparse
import logging
import os
import sys
import pyds9
from astropy import units
from astropy.time import Time
from ossos import storage
from ossos.wcs import WCS
import tempfile

from . import settings

config = settings.AppConfig()


def start_ds9(name):
    """
    Start a DS9 program with target 'name'
    :param name: name to assign to DS9 as target for xpaset commands
    :return:
    """

    ds9 = pyds9.DS9(target=name)
    levels = ['INIT', 'PREF']
    for level in levels:
        setting = config.read(f"DS9.{level}")
        for key in list(setting.keys()):
            ds9.set(f"{key.replace('_', ' ')} {setting[key]}")
    ds9.set("frame delete all")
    return ds9


def get_ds9(name):
    """
    Get a connection to DS9 via the name.

    :param name:
    :return:
    """
    return pyds9.DS9(target=name, start=False)


def main(**kwargs):
    """
    load images from vos:OSSOS/dbimages as referenced from .ast lines.

    :param kwargs:

    Exect to the 'orbit' and 'ossos_observations'
    :return:
    """
    from daomop import daophot
    orbit = kwargs['orbit']
    observations = kwargs['ossos_observations']

    # orbit = kwargs.get('orbit', None)
    # isinstance(BKOrbit, orbit)

    ds9 = get_ds9('validate')
    # Load the 3 images associated with this point/ccd/rate/angle set.

    wcs_dict = {}
    ds9.set('frame delete all')
    ds9.set('zscale')
    displayed_images = []
    for frame in observations:
        observation = observations[frame]
        observation: Observation
        expnum, ccd = frame.split('p')
        hdulist = None
        while True:
            try:
                hdulist = storage.ra_dec_cutout(storage.get_uri(expnum),
                                                observation.coordinate,
                                                radius=10 * units.arcsec)
            except Exception as ex:
                logging.warning(ex)
                continue
            break
        if hdulist is None:
            logging.error(f"Failed to get image for {frame} at {observation.coordinate}")
            continue
        # Get the extension that has the RA/DEC location
        # and also delete all the PV keywords once we have WCS object.
        ra = observation.coordinate.ra.degree
        dec = observation.coordinate.dec.degree
        x = y = None
        for hdu in hdulist:
            if hdu.header['NAXIS'] != 2:
                continue
            wcs_dict[frame] = WCS(hdu.header)
            del(hdu.header['PV*'])
            # skip primary header
            x, y = wcs_dict[frame].sky2xy(ra, dec)
            if 0 < x < hdu.header['NAXIS1'] and 0 < y < hdu.header['NAXIS2']:
                break
        if x is None or y is None:
            logging.error(f"Failed to get x/y coordinate from ra/dec {ra:9.5},{dec:9.5} for frame {frame}")
            continue
        with tempfile.NamedTemporaryFile() as fobj:
            hdulist.writeto(fobj.name, overwrite=True)
            mjdmid = (hdulist[1].header['MJD-OBS'] + hdulist[1].header['MJDEND']) / 2.0
            obsdate = Time(mjdmid, format='mjd')
            orbit.predict(obsdate)
            mark_ra, mark_dec = wcs_dict[frame].xy2sky(x, y, usepv=False)
            uncertainty_ellipse = (orbit.dra.to('arcsec').value,
                                   orbit.ddec.to('arcsec').value,
                                   orbit.pa.to('degree').value + 90)
            ds9.set('frame new')
            ds9.set(f'mosaicimage wcs {fobj.name}')
            # ds9.set_pyfits(hdulist)
            displayed_images.append(frame)
            ds9.set('regions', f'icrs; ellipse({mark_ra.value},{mark_dec.value},'
                               f'{uncertainty_ellipse[0]*5}",'
                               f'{uncertainty_ellipse[1]*5}",'
                               f'{uncertainty_ellipse[2]})')
            ds9.set(f'pan to {ra} {dec} wcs icrs')

    if not len(ds9.get('frame')) > 0:
        return {}

    ds9.set('frame first')
    ds9.set('frame center')
    ds9.set('frame match image')
    obs = {}
    images = displayed_images
    # Build a map of allowed key strokes
    allowed_keys = {'x': ('x', 'centroid at this location'),
                    'q': ('q', 'Quit this image set'),
                    'Q': ('Q', 'Quit the session'),
                    'p': ('p', 'Previous frame'),
                    'n': ('n', 'Next Frame'),
                    'r': ('r', 'Create a NULL observation'),
                    'C': ('C', 'Centroid frame and match wcs')}

    for key in [x.split() for x in config.read("MPC.NOTE1OPTIONS")]:
        allowed_keys[key[0]] = key

    while True:
        try:
            result = ds9.get('iexam key coordinate wcs icrs degrees')
            key, mark_ra, mark_dec = result.split()
            mark_ra = float(mark_ra)
            mark_dec = float(mark_dec)
            logging.debug(f"DS9 Returned: {result} -> {key} {mark_ra} {mark_dec}")
        except Exception as ex:
            logging.debug(f"DS9 get exception: {ex}")
            continue

        if key not in allowed_keys:
            logging.info(f"Allowed keys: ")
            for key in allowed_keys:
                print(f"{key} -> {allowed_keys[key][1]}")
            continue

        if key == 'n':
            ds9.set('frame next')
            continue
        if key == 'p':
            ds9.set('frame prev')
            continue
        if key == 'C':
            ds9.set('frame center')
            ds9.set('frame match wcs')
            continue
        if key == 'q':
            break
        if key == 'Q':
            sys.exit(0)

        note1 = allowed_keys[key][0]
        frame_no = int(ds9.get('frame')) - 1
        frame = images[frame_no]
        x, y = wcs_dict[frame].sky2xy(mark_ra, mark_dec, usepv=False)
        ds9.set('regions', f'image; circle {x} {y} 20')
        hdulist = ds9.get_pyfits()
        hdulist[0].writeto('temp.fits', overwrite=True)
        centroid = not note1 == 'H'
        phot = daophot.phot_mag('temp.fits',
                                [x, ], [y, ],
                                aperture=5,
                                sky_inner_radius=15,
                                sky_annulus_width=10,
                                apcor=0.3,
                                zmag=27.1,
                                maxcount=1000,
                                extno=0,
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

        mjdmid = (hdulist[0].header['MJD-OBS'] + hdulist[0].header['MJDEND']) / 2.0
        obsdate = Time(mjdmid, format='mjd', precision=6).mpc

        try:
            ra, dec = wcs_dict[frame].xy2sky(cen_x, cen_y)
            ds9.set('regions', f'icrs; circle({ra},{dec},0.2")')
            logging.debug(f"Got {ra},{dec} from {cen_x},{cen_y}")
        except Exception as ex:
            logging.warning(f"Failure converting {cen_x, cen_y} to RA/DEC for {frame}")
            logging.warning(ex)
            continue

        record_key = os.path.basename(frame)
        observation = observations[frame]
        obs[record_key] = (Observation(
            null_observation=key == 'r',
            provisional_name=observation.provisional_name,
            note1=note1,
            note2='C',
            date=obsdate,
            ra=ra * units.degree,
            dec=dec * units.degree,
            mag=obs_mag,
            mag_err=obs_mag_err,
            band='r',
            observatory_code='568',
            comment=None,
            xpos=x,
            ypos=y,
            frame=os.path.basename(frame),
            astrometric_level=2))
        ds9.set('frame next')
    return obs


def _main(**kwargs):
    """
    Start up ds9, load the ephemeris file and then measure the sources.

    :param kwargs:
    :return:
    """
    start_ds9('validate')
    ast_filename = kwargs['ast_filename']

    logging.info(f"Attempting measures of {ast_filename}")
    all_observations = EphemerisReader().read(ast_filename)
    orbit = BKOrbit(all_observations)
    logging.info(orbit.summarize())

    # Select those observations with FRAME information
    ossos_observations = {}
    non_ossos_observations = []
    observation: Observation
    for observation in all_observations:
        try:
            record_key = observation.comment.frame
            if len(record_key) > 0:
                if record_key in ossos_observations:
                    logging.error(f"Duplicate frame value on two or more ast lines: {record_key}. Keeping first line.")
                    continue
                ossos_observations[record_key] = observation
            else:
                non_ossos_observations.append(observation)
        except Exception as ex:
            logging.warning(f"{ex}: {observation}")
            continue

    logging.info(f"Inspecting {len(orbit.observations)} images.")
    kwargs['orbit'] = orbit
    kwargs['ossos_observations'] = ossos_observations
    new_obs = main(**kwargs)

    # replace existing measurement lines with new ones.
    logging.info(f"{new_obs}")
    for record_index in new_obs:
        ossos_observations[record_index] = new_obs[record_index]

    # save beside input .ast file
    output_ast_filename = ast_filename + "_vetted"
    with open(output_ast_filename, 'w') as mpc_obj:
        for record in ossos_observations:
            mpc_obj.write(ossos_observations[record].to_string() + "\n")
        for observation in non_ossos_observations:
            mpc_obj.write(observation.to_string() + "\n")
    try:
        orbit = BKOrbit(None, ast_filename=output_ast_filename)
        logging.info(orbit.summarize())
    except Exception as ex:
        logging.error(f"{ex}")


def get_valid_obs_count(observations):
    """
    Count the number of observation that are not null observations.
    :param observations:
    :return:
    """
    nobs = 0
    for obs in observations:
        if not obs.null_observation:
            nobs += 1
    return nobs


def run():
    """
    Hook to allow make this into a consol script.
    :return:
    """
    main_parser = argparse.ArgumentParser()
    main_parser.add_argument('ast_filename', type=str)
    main_parser.add_argument('--log-level', choices=['DEBUG', 'INFO', 'ERROR'], default='INFO')
    main_parser.add_argument('--skip', action='store_true')
    args = main_parser.parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level))
    _main(ast_filename=args.ast_filename, skip=args.skip)


if __name__ == '__main__':
    run()
