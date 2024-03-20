"""
Create an Observation using ds9 displaying an image of a KBO source.
"""
import argparse
import logging
import os, sys
from mp_ephem import BKOrbit
import numpy
import math
import pyds9
from astropy import units
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from astropy.wcs import WCS
from mp_ephem import BKOrbit, EphemerisReader
from mp_ephem.ephem import Observation
from vos import Client
from tempfile import NamedTemporaryFile

from . import settings
from . import util

config = settings.AppConfig()



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


def main(**kwargs):
    """

    :param kwargs:
    :type orbit: BKOrbit
    :return:
    """
    from daomop import daophot
    orbit = kwargs['orbit']
    images = kwargs['images']

    # orbit = kwargs.get('orbit', None)
    # isinstance(BKOrbit, orbit)

    ds9 = get_ds9('validate')
    # Load the 3 images associated with this point/ccd/rate/angle set.

    wcs_dict = {}
    offset = {}
    ds9.set('frame delete all')
    ds9.set('zscale')
    displayed_images = []
    for image in images:
        with fits.open(image) as hdulist:
            header = hdulist[1].header
            obsdate = Time(hdulist[0].header['DATE-AVG'], scale='tai')
            try:
                wcs_header_filename = image.replace('.fits','.mega.head')
                wcs_header = fits.Header.fromtextfile(wcs_header_filename)
                wcs_dict[image] = WCS(wcs_header)
            except Exception as ex:
                wcs_header = header
                wcs_header_filename = image
                wcs_dict[image] = WCS(wcs_header)
                logging.debug(f"using original wcs")
            if orbit is not None:
                orbit.predict(obsdate)
                ra = orbit.coordinate.ra.degree
                dec = orbit.coordinate.dec.degree
                rad = int(max(orbit.dra.to('arcsec').value, orbit.ddec.to('arcsec').value)/0.17)
                uncertainty_ellipse = (orbit.dra.to('arcsec').value,
                                       orbit.ddec.to('arcsec').value,
                                       orbit.pa.to('degree').value + 90)
            else:
                uncertainty_ellipse = 3, 3, 0
                rad = int(3/0.17)
            x, y = wcs_dict[image].all_world2pix(ra, dec, 0)
            if x < 30 or x > 2048 - 30 or y < 30 or y > 4176 - 30:
                logging.warning(f"Skipping (image) as too near chip edge")
                continue
            cutsize = max(100, 3*rad)
            x1 = int(max(0, x-cutsize))
            x2 = int(min(header['NAXIS1'], x+cutsize))
            y1 = int(max(0, y-cutsize))
            y2 = int(min(header['NAXIS2'], y+cutsize))
            offset[image] = x1, y1 
            wcs_header['CRPIX1'] -= offset[image][0]
            wcs_header['CRPIX2'] -= offset[image][1]
            display_hdu = fits.HDUList([fits.PrimaryHDU(data=hdulist[1].data[y1:y2,x1:x2],
                                                        header=wcs_header)])
            ds9.set('frame new')
            displayed_images.append(image)
            ds9.set_pyfits(display_hdu)
            ds9.set('regions', f'icrs; ellipse({ra},{dec},'
                               f'{uncertainty_ellipse[0]}",'
                               f'{uncertainty_ellipse[1]}",'
                               f'{uncertainty_ellipse[2]})')
            ds9.set(f'pan to {ra} {dec} wcs icrs')
    if not len(ds9.get('frame'))> 0:
        return {}
        
    ds9.set('frame match wcs')
    ds9.set('frame first')
    obs = {}
    images = displayed_images
    # Build a map of allowed key strokes
    allowed_keys = {'x': ('', 'centroid at this location'),
                    'q': ('', 'Quit this image set'),
                    'Q': ('', 'Exit the program'),
                    'p': ('', 'Previous frame'),
                    'n': ('', 'Next Frame'),
                    'r': ('', 'Create a NULL observation')}
    
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
            ds9.set('frame prev')
            continue

        if key not in allowed_keys:
            logging.info(f"Allowed keys: ")
            for key in allowed_keys:
                print(f"{key} -> {allowed_keys[key][1]}")
            continue

        if key == 'q':
            break
        if key == 'Q':
            sys.exit(0)
        note1 = allowed_keys[key][0]
        frame_no = int(ds9.get('frame')) - 1
        image = images[frame_no]
        with fits.open(image) as _hdulist:
            exptime = _hdulist[0].header.get('EXPTIME', 0.0)
        ds9.set('regions', f'image; circle {x} {y} 20 # color=blue ')
        hdulist = ds9.get_pyfits()
        with NamedTemporaryFile(mode='w+b', delete=False, suffix=".fits") as fobj:
            hdulist[0].writeto(fobj.name, overwrite=True)
            centroid = not note1 == 'H'
            phot = daophot.phot_mag(fobj.name,
                                    [x, ], [y, ],
                                    aperture=5,
                                    sky_inner_radius=15,
                                    sky_annulus_width=10,
                                    apcor=0.3,
                                    zmag=26.7,
                                    maxcount=1000,
                                    extno=0,
                                    exptime=exptime,
                                    centroid=centroid)
            os.unlink(fobj.name)
        
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

        with fits.open(image) as _hdu_list:
            primary_hdu = _hdu_list[0]
            obsdate = Time(Time(primary_hdu.header['DATE-AVG'], scale='tai').mjd,
                       format='mjd',
                       precision=6).mpc
            frame = primary_hdu.header.get('FRAMEID', os.path.basename(image))
        try:
            ra, dec = wcs_dict[image].all_pix2world(cen_x + offset[image][0],
                                                    cen_y + offset[image][1],
                                                    1)
            ds9.set('regions', f'icrs; circle({ra},{dec},0.2")')
            logging.debug(f"Got {ra},{dec} from {cen_x},{cen_y}")
        except Exception as ex:
            logging.warning(f"Failure converting {cen_x, cen_y} to RA/DEC for {image}")
            logging.warning(ex)
            logging.warning(f"Got: {ra},{dec}")

        # record_key = os.path.basename(image)
        record_key = obsdate
        obs[record_key] = (Observation(
            null_observation=key == 'r',
            provisional_name=kwargs['provisional_name'],
            note1=note1,
            note2='C',
            date=obsdate,
            ra=ra*units.degree,
            dec=dec*units.degree,
            mag=obs_mag,
            mag_err=None,
            band='r',
            observatory_code='568',
            comment=None,
            xpos=x,
            ypos=y,
            frame=frame,
            astrometric_level=2))
        discovery = False
        ds9.set('frame next')
    return obs


def _main(**kwargs):
    start_ds9('validate')
    ast_filename = kwargs['ast_filename']
    output_ast_filename = ast_filename
    
    logging.info(f"Attempting measures of {kwargs['provisional_name']}, will write to {ast_filename}")
    obs = {}
    kwargs['orbit'] = BKOrbit(None, ast_filename)
    unique_obs = []
    for ob in kwargs['orbit'].observations:
        record_key = ob.date.mpc
        #try:
        #    record_key = ob.comment.frame
        #    if len(record_key) > 0:
        #        if record_key in obs:
        #            logging.warning(f"Duplicate frame value: {record_key}")
        #            continue
        #except:
        #    pass
        obs[record_key] = ob

    orb = BKOrbit([ obs[x] for x in obs])
    logging.info(orb.summarize())
    logging.info(f"Measuring on {len(kwargs['images'])} images, {kwargs['nframes']} at a time.")
    step_size = kwargs['nframes']
    stride = kwargs['stride']
    niters = int(math.ceil(len(kwargs['images']) / (stride*step_size) ))
    images = kwargs['images']
    for i in range(niters):
        image_set = images[i*step_size*stride:min(len(kwargs['images']),i*step_size*stride+step_size*stride):stride]
        if kwargs['skip']:
            kwargs['images'] = []
            for image in image_set:
                frame = fits.open(image)[0].header.get('FRAMEID', os.path.basename(image)[0:12])
                if frame in obs:
                    continue
                kwargs['images'].append(image)
        else:
            kwargs['images'] = image_set
            
        if not len(kwargs['images']) > 0:
            continue
        kwargs['orbit'] = BKOrbit(None, ast_filename=output_ast_filename)
        new_obs = main(**kwargs)
        logging.debug(f"{new_obs}")
        for record_index in new_obs:
            obs[record_index] = new_obs[record_index]

        with open(output_ast_filename, 'w') as mpc_obj:
            for record in obs:
                mpc_obj.write(obs[record].to_string() + "\n")
        try:
            orbit = BKOrbit(None, ast_filename=output_ast_filename)
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
    kwargs['p_name'] = util.get_provisional_name(**kwargs)
    _main(**kwargs)


def run():
    main_parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    main_parser.add_argument('provisonal_name')
    main_parser.add_argument('ast_filename', type=str)
    main_parser.add_argument('images', nargs='+')
    main_parser.add_argument('--nframes', type=int, default=10, help="Number of frames from images list to do in a sequence")
    main_parser.add_argument('--stride', type=int, default=1, help="Skip this number of images per step while doing list of images")
    main_parser.add_argument('--log-level', choices=['DEBUG', 'INFO', 'ERROR'], default='INFO')
    main_parser.add_argument('--skip', action='store_true', help="Skip images whose frame_id values are already in the astrometry input file.")
    args = main_parser.parse_args()
    _format="%(asctime)s :: %(levelname)s :: %(module)s.%(funcName)s:%(lineno)d %(message)s"
    if args.log_level == 'INFO':
         _format="%(message)s"
    logging.basicConfig(level=getattr(logging, args.log_level), format=_format)
    _main(images=args.images, ast_filename=args.ast_filename, provisional_name=args.provisonal_name, nframes=args.nframes, stride=args.stride, skip=args.skip)


if __name__ == '__main__':
    run()
