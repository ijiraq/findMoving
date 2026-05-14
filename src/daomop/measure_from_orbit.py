"""
Create an Observation using ds9 displaying an image of a KBO source.
"""
import argparse
import sys
import os
import logging
import math
from typing import List

import pyds9
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from astropy.wcs import WCS
from astropy import units
from mp_ephem import BKOrbit, EphemerisReader
from . import settings
from . import util
from . import daophot
from .fwhm import fit_fwhm, aperture_correction
from mp_ephem.ephem import Observation
from dataclasses import dataclass, field

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


def mark_planted_sources(ds9, header):

    try:
        plant_list_column_names = "index ra dec x y rate angle rate_ra rate_dec mag psf_amp g_i".split()
        WARPDIR = header['WARPD001']
        PLANTFILE = header['PLANT001']
        plant_table = Table.read(f"{WARPDIR}/{PLANTFILE}", format='ascii', names=plant_list_column_names)
        for row in plant_table:
            ra = row['ra']
            dec = row['dec']
            ds9.set('regions', f'icrs; circle({ra},{dec},0.2") # color=cyan width=2')
    except Exception as ex:
        logging.warning(f"Failed to mark planted sources: {ex}")
    return


def create_observation_record(image:str, provisional_name:str, x:float, y:float, note1:str,
                              aperture:float, apcor:float, comment:str='stack', likelihood:int=-1,
                              null_observation:bool=False) -> Observation:
    """
    Create an Observation record from the given hdulist, RA/DEC, and other context.
    """
    obs_mag = None
    obs_mag_err = None
    cen_x = x
    cen_y = y
    filter_band_map = {'r2': 'r', 'default': 'r'}
    observatory_codes = {'CFHT': 568, 'default': 568}
    if not null_observation:
        try:
            centroid = not note1.lower() == 'h'
            phot = daophot.phot_mag(image,
                                    [x, ], [y, ],
                                    aperture=aperture,
                                    sky_inner_radius=15,
                                    sky_annulus_width=10,
                                    apcor=apcor,
                                    zmag=None,
                                    maxcount=10000,
                                    extno=1,
                                    centroid=centroid)
            phot_failure = (phot['PIER'][0] != 0 or phot.mask[0]['MAG'] or phot.mask[0]['MERR'])
            sky_failure = phot['SIER'][0] != 0
            cen_failure = phot['CIER'][0] != 0
            if phot_failure or sky_failure or cen_failure:
                logging.warning(f"iraf.daophot.phot error:\n {phot}")
                note1 = "H"
            else:
                cen_x = phot['XCENTER'][0]
                cen_y = phot['YCENTER'][0]
                obs_mag = phot['MAG'][0]
                obs_mag_err = phot['MERR'][0]
        except Exception as ex:
            logging.warning(f"Photometry failed: {ex}")
            note1 = "H"

    with fits.open(image) as hdulist:
        wcs = WCS(hdulist[1].header)
        band = hdulist[0].header.get('FILTER', 'default')
        band = filter_band_map.get(band, filter_band_map['default'])
        observatory_code = observatory_codes.get(hdulist[0].header.get('ORIGIN', 'CFHT'),
                                                 observatory_codes['default'])
        astrometric_level = hdulist[0].header.get('ASTLEVEL', 0)
        xoffset = hdulist[0].header.get('XOFFSET', 0.0)
        yoffset = hdulist[0].header.get('YOFFSET', 0.0)
        obsdate = Time(Time(hdulist[0].header['DATE-AVG'], scale='tai').mjd, format='mjd', precision=5).mpc
        frame_val = hdulist[0].header.get('FRAMEID', os.path.basename(image))
        ra_val, dec_val = wcs.all_pix2world(cen_x, cen_y, 1)
    return Observation(
        discovery=False,
        likelihood=likelihood,
        survey_code='C',
        null_observation=null_observation,
        provisional_name=provisional_name,
        note1=note1,
        note2='C',
        date=obsdate,
        ra=ra_val*units.degree,
        dec=dec_val*units.degree,
        mag=obs_mag,
        mag_err=obs_mag_err,
        band=band,
        observatory_code=observatory_code,
        comment=f"{comment} {aperture:.2f}:{apcor:.2f}",
        xpos=cen_x+xoffset,
        ypos=cen_y+yoffset,
        frame=frame_val,
        astrometric_level=astrometric_level
    )

@dataclass
class FakeDS9:
    """
    A fake DS9 class to use when we don't have a ds9 instance.
    """
    _x: List[float] = field(default_factory=list)
    _y: List[float] = field(default_factory=list)
    _frame: int = -1

    @property
    def frame(self) -> int:
        return self._frame

    @frame.setter
    def frame(self, value: str|int):
        if "delete all" in value.lower():
            self._x.clear()
            self._y.clear()
            self._frame = -1
            return
        steps = {'next': self._frame+1,
                 'prev': self._frame-1,
                 'new': self._frame+1,
                 'first': 0,
                 'last': len(self._x)-1}
        try:
            self._frame = int(value)
        except ValueError:
            self._frame = steps.get(value.lower(),
                                    self._frame)

    @property
    def x(self) -> float:
        if self._frame < 0 or self._frame >= len(self._x):
            return None
        return self._x[self._frame]

    @property
    def y(self) -> float:
        if self._frame < 0 or self._frame >= len(self._y):
            return None
        return self._y[self._frame]

    @x.setter
    def x(self, value: str|float|int):
        """
        Set the x coordinate for the current frame.
        :param value: The x coordinate as a string.
        """
        if self._frame < 0:
            logging.warning("FakeDS9.set_x called with no frame set.")
            return
        try:
            self._x[self._frame] = float(value)
        except IndexError:
            self._x.append(float(value))
        except ValueError:
            logging.error(f"Invalid x value: {value}")

    @y.setter
    def y(self, value: str|float|int):
        """
        Set the y coordinate for the current frame.
        :param value: The y coordinate as a string.
        """
        if self._frame < 0:
            logging.warning("FakeDS9.set_y called with no frame set.")
            return
        try:
            self._y[self._frame] = float(value)
        except IndexError:
            self._y.append(float(value))
        except ValueError:
            logging.error(f"Invalid y value: {value}")

    def set(self, *args, **kwargs):
        values = args[0].split()
        if not hasattr(self, values[0]):
            logging.debug(f"FakeDS9.set called with unknown attribute: {values[0]}")
            return
        setattr(self, values[0], values[1])  # Set the attribute based on the first value

    def set_pyfits(self, *args, **kwargs):
        logging.debug(f"FakeDS9.set_pyfits called with args: {args}, kwargs: {kwargs}")

    def get(self, *args, **kwargs):
        logging.debug(f"FakeDS9.get called with args: {args}, kwargs: {kwargs}")
        values = args[0].split()
        if 'iexam' in values[0]:
            if self._frame > len(self._x) - 1:
                return f'q 0.0 0.0' # No valid frame, return a default value
            return f'a {self.x} {self.y}'  # Simulate a key press at (x, y)
        if 'frame' in values[0]:
            return self._frame + 1 # Return the current frame number (1-indexed)
        return None

def main(**kwargs):
    """

    :param kwargs:
    :type orbit: BKOrbit
    :return:
    """
    orbit = kwargs['orbit']
    images = kwargs['images']
    aperture = kwargs['aperture']
    apcor = kwargs['apcor']
    auto_process = kwargs.get('auto_process', False) # if True do not display the images, just process them.

    # orbit = kwargs.get('orbit', None)
    # isinstance(BKOrbit, orbit)

    if not auto_process:
        ds9 = get_ds9('validate')
    else:
        ds9 = FakeDS9()
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
                logging.debug(f"using wcs in {wcs_header_filename}")
            except Exception as ex:
                wcs_header = header
                wcs_dict[image] = WCS(wcs_header)
                logging.debug(f"using wcs in {image}")
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
            cutsize = max(100, 3*rad)
            x, y = wcs_dict[image].all_world2pix(ra, dec, 0)
            if x < -cutsize or x > 2048+cutsize or y < -cutsize or y > 4176+cutsize :
                logging.warning(f"Skipping {image}: predicted source location ({x},{y}) +/- ({rad}) off image")
                continue
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
            ds9.set('contour smooth 6')
            ds9.set('contour nlevels 7')
            ds9.set('contour generate')
            ds9.set('contour yes')
            if auto_process:
                ds9.set(f'x {x-offset[image][0]}')
                ds9.set(f'y {y-offset[image][1]}')
            ds9.set('regions', f'icrs; ellipse({ra},{dec},'
                               f'{uncertainty_ellipse[0]}",'
                               f'{uncertainty_ellipse[1]}",'
                               f'{uncertainty_ellipse[2]}) # color=red width=2')
            mark_planted_sources(ds9, hdulist[0].header)
            ds9.set(f'pan to {ra} {dec} wcs icrs')

    try:
        frameno = int(ds9.get('frame'))
    except Exception as ex:
        logging.error(f"Failed to get frame number from DS9: {ex}")
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
        null_observation = key.lower() == 'r'
        note1 = allowed_keys[key][0]
        frame_no = int(ds9.get('frame')) - 1
        image = images[frame_no]
        # put x,y from the ds9 coordinate system into the image file coordinate system.
        x += offset[image][0]
        y += offset[image][1]
        if kwargs['auto_apcor']:
            fwhm = 5
            with fits.open(image) as hdulist:
                try:
                    fwhm = fit_fwhm(hdulist[1].data, (x-offset[image][0], y-offset[image][1]))
                except ValueError as ve:
                    logging.warning(f"Failed to measure FWHM for {image} at ({x}, {y}): {ve} using default {fwhm}")
            aperture = 1.1 * fwhm
            apcor = aperture_correction(fwhm, aperture)
        obs_record = create_observation_record(image,
                                               kwargs['provisional_name'], x, y, note1,
                                               aperture, apcor, comment='stack', null_observation=null_observation)
        ds9.set('regions', f'image; circle {x-offset[image][0]} {y-offset[image][1]} 4 # color=blue width=2')
        record_key = obs_record.date.mpc[0:15]
        obs[record_key] = obs_record
        ds9.set('frame next')
        update_ast_file(obs, kwargs['tlf_filename'])
    return obs


def update_ast_file(obs: dict, output_ast_filename: str) -> None:
    old_obs = {}
    if os.access(output_ast_filename, os.F_OK):
        for ob in  EphemerisReader().read(output_ast_filename):
            old_obs[ob.date.mpc[0:15]] = ob
    old_obs.update(obs)
    with open(output_ast_filename, 'w') as mpc_obj:
        for record in old_obs:
            mpc_obj.write(old_obs[record].to_tnodb() + "\n")
    return


def _main(**kwargs):
    if not kwargs.get('auto_process', False):
        start_ds9('validate')
    ast_filename = kwargs['ast_filename']
    output_ast_filename = ast_filename
    
    logging.info(f"Attempting measures of {kwargs['provisional_name']}, will write to {ast_filename}")
    kwargs['orbit'] = BKOrbit(None, ast_filename)
    kwargs['auto_apcor'] = kwargs.get('auto_apcor', False)
    obs = {}
    for ob in kwargs['orbit'].observations:
        obs[ob.date.mpc] = ob

    orb = BKOrbit([obs[x] for x in obs])
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
        kwargs['tlf_filename'] = f"{kwargs['provisional_name']}.inp"
        new_obs = main(**kwargs)
        logging.debug(f"{new_obs}")
        # Don't overwrite the previous astrometry in this code. 
        # for record_index in new_obs:
        #     obs[record_index] = new_obs[record_index]
        # update_ast_file(obs, output_ast_filename)


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
    main_parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                          fromfile_prefix_chars='@')
    main_parser.add_argument('provisonal_name')
    main_parser.add_argument('ast_filename', type=str)
    main_parser.add_argument('images', nargs='+')
    main_parser.add_argument('--nframes', type=int, default=10, help="Number of frames from images list to do in a sequence")
    main_parser.add_argument('--stride', type=int, default=1, help="Skip this number of images per step while doing list of images")
    main_parser.add_argument('--log-level', choices=['DEBUG', 'INFO', 'ERROR'], default='INFO')
    main_parser.add_argument('--skip', action='store_true', help="Skip images whose frame_id values are already in the astrometry input file.")
    main_parser.add_argument('--apcor', type=float, help="Set the aperture correction value", default=0.25)
    main_parser.add_argument("--auto-apcor", action='store_true',
                             help="Automatically determine the aperture size aperture correction from fwhm")
    main_parser.add_argument('--aperture', type=float, help="Aperture to measure flux with", default=5.0)
    main_parser.add_argument('--photzp', type=str, help="Header keyword with zeropoint", default='PHOTZP')
    main_parser.add_argument('--auto-process', action='store_true', help="Do not display images, just process them.")

    args = main_parser.parse_args()
    _format="%(asctime)s :: %(levelname)s :: %(module)s.%(funcName)s:%(lineno)d %(message)s"
    if args.log_level == 'INFO':
         _format="%(message)s"
    logging.basicConfig(level=getattr(logging, args.log_level), format=_format)
    _main(images=args.images, ast_filename=args.ast_filename, provisional_name=args.provisonal_name,
          nframes=args.nframes, stride=args.stride, skip=args.skip, apcor=args.apcor,
          aperture=args.aperture, photzp=args.photzp, auto_process=args.auto_process, autoapcor=args.auto_apcor)


if __name__ == '__main__':
    run()
