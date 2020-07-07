import argparse
import glob
import logging
import os

import numpy as np
from astropy import time, units
from astropy.io import fits
from reproject import reproject_interp
from reproject.mosaicking import find_optimal_celestial_wcs, reproject_and_coadd


STACKING_MODES = {'MEDIAN': np.nanmedian,
                  'MEAN': np.nanmean,
                  'SUM': np.nansum,
                  'MAX': np.nanmax,
                  'DEFAULT': np.nanmedian}

def swarp(hdu, wcs_out, shape_out):
    """
    Before we can stack the data we need to project onto an aligned grid.
    """
    reproject = reproject_interp
    new_header = hdu.header.copy()
    wcs = {'CTYPE1': 'RA---TAN', 'CTYPE2': 'DEC--TAN',
           'CD1_1': pixel_scale, 'CD1_2': 0.0, 'CD2_2': pixel_scale, 'CD2_1': 0.0}
    for key in wcs:
        new_header[key]=wcs[key]
    return reproject(hdu, new_header, return_footprint=False)

    

def shift(data, mjd, rx, ry, pix_scale, rf=3, stacking_mode='MEDIAN', mid_mjd=None):
    """

    :param data: array of fits data
    :type data: np.array
    :param mjd: array of observation centre times (same order as datas)
    :type mjd: np.array
    :param rx: rate in 'dx/dt' to shift data array at in units of arcsec/hour
    :type rx: float
    :param ry: date in 'dy/dy' to shift data array at in units of arcsec/hour
    :type ry: float
    :param pix_scale: scale of image pixels in arcsec/pixel
    :type pix_scale: float
    :param rf: rescaling factor, up-sample by this factor before shifing at integer amounts.
    :type rf: int
    :param stacking_mode: algorithm used to combine pixel values, default: MEDIAN (MEAN, MEADIN, SUM)
    :type stacking_mode: str
    :param mid_mjd:  MJD time at centre of stack (taken as the mean of mjds if None)
    :type mid_mjd: float
    :rtype: np.array
    :return: combined datas after shifting at dx/dy and combined using stacking_mode.
    """
    from trippy.trippy_utils import downSample2d

    stacking_mode = STACKING_MODES.get(stacking_mode, STACKING_MODES['DEFAULT'])

    print('Shifting at {} {} "/hr'.format(rx, ry))

    if mid_mjd is None:
        mid_mjd = np.mean(mjd)

    # which of the mjd in the array is closest to the mid_mjd value.
    mid = np.argmin(np.abs(mjd - mid_mjd))

    # Hold that each up-sampled data array will have
    (A, B) = data[0].shape
    A *= rf
    B *= rf

    # time since mid time.
    dt = (mjd - mid_mjd) * 24.0
    # amount to shift, scaled for up-sampling and converted to integer pixels from arcsec/hour
    dx = (-rx * dt / pix_scale * rf + 0.5).astype('int')
    dy = (-ry * dt / pix_scale * rf + 0.5).astype('int')

    logging.debug("Up-sampled x-pixel shift is:{}".format(dx))
    logging.debug("Up-sampled y-pixel shift rate is:{}".format(dy))

    # outs contains the shifted versions of the arrays after down sampling.
    outs = []
    for ii in range(len(data)):
        logging.debug("Working on array {} of size {} and shifting by dx {} and dy {}".format(ii + 1,
                                                                                              len(data),
                                                                                              dx[ii],
                                                                                              dy[ii]))
        # If this is the array that is the mid point of the stack
        # then just add to the data stack as no shifting is needed.
        if ii == mid or (dy[ii] == 0 and dx[ii] == 0):
            outs.append(data[mid])
            continue

        # up-scale the array onto an array that is 'rf' times bigger.
        rep = np.repeat(np.repeat(data[ii], rf, axis=0), rf, axis=1)
        logging.debug("Data from shape {} has been sampled into shape {}".format(data.shape,
                                                                                 rep.shape))
        # a/b (c/d) are the range of the first (second) index where the data should go into
        # aa/bb (cc/dd) are the range of the first (second) index where the data come from.
        if dy[ii] < 0:
            a, b = 0, A + dy[ii]
            aa, bb = -dy[ii], A
        elif dy[ii] > 0:
            a, b = dy[ii], A
            aa, bb = 0, A - dy[ii]
        else:
            a, b = 0, A
            aa, bb = 0, A

        if dx[ii] < 0:
            c, d = 0, B + dx[ii]
            cc, dd = -dx[ii], B
        elif dx[ii] > 0:
            c, d = dx[ii], B
            cc, dd = 0, B - dx[ii]
        else:
            c, d = 0, B
            cc, dd = 0, B

        rep[a:b, c:d] = rep[aa:bb, cc:dd]
        outs.append(downSample2d(rep, rf))

    return stacking_mode(outs, axis=0)

def shift_rates(r_min, r_max, angle_min, angle_max):
    """
    @param rmin: minimum shift rate (''/hour)
    @param rmax: maximum shift rate (''/hour)
    @param angle_min: minimum angle to shift at (degrees)
    @param angle_max: maximum angle to shift at (degrees)
    """
    rates = []
    for dr in np.linspace(r_min, r_max, int((r_max - r_min) / 0.25) + 1):
        for dd in np.linspace(angle_min, angle_max, int((angle_max-angle_min) / .25) + 1):
            rates.append([dr*units.arcsecond/units.hour, dd*units.degree])
            logging.debug(f'Rate: {dr} and Angel:{dd}')
    return rates

def mid_exposure_mjd(hdu):
    
    mjd_start = time.Time(hdu.header['MJD-STR'], format='mjd')
    mjd_end = time.Time(hdu.header['MJD-END'], format='mjd')
    return mjd_start + (mjd_end - mjd_start)/2.0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('basedir', help="Root directory of LSST pipelined data")
    parser.add_argument('--pointing', help="sky patch to process (e.g. 0,0)", nargs=1)
    parser.add_argument('--rerun', help="rerun directory containing the warped difference images.", nargs=1)
    parser.add_argument('--filter', help="Filter to stack", default="HSC-R2")
    parser.add_argument('--pixel-scale', help="What should the pixel scale of the stack be? (in arcsecond)", 
                        default=0.16)
    parser.add_argument('--ccd', help="Which CCD to stack?", type=int, default=0)
    parser.add_argument('--exptype', help="Whay type of exposures to coadd?", default='deepDiff')
    parser.add_argument('--log-level', help="What level to log at? (ERROR, INFO, DEBUG)", default="ERROR",
                        choices=['INFO', 'ERROR', 'DEBUG'])

    args = parser.parse_args()
    levels = {'INFO': logging.INFO, 'ERROR': logging.ERROR, 'DEBUG': logging.DEBUG}
    logging.basicConfig(level=levels[args.log_level])
    pixel_scale = args.pixel_scale*units.arcsec
    ccd = f'{args.ccd:03d}'
    reruns = args.rerun[0].split(":")
    if len(reruns) > 2:
        raise ValueError("Don't know what to do with more then 2 rerun directories.")

    if len(reruns) == 1:
        input_rerun = output_rerun = reruns[0]
    else:
        input_rerun = reruns[0]
        output_rerun = reruns[1]

    input_dir = os.path.join(args.basedir, 'rerun', input_rerun, args.exptype, 
                             args.pointing[0], args.filter, f'DIFF*-{ccd}.fits')
    logging.info(f'Loading all images matching pattern: {input_dir}')
    images = glob.glob(input_dir)
    if not len(images) > 0:
        raise OSError(f'No images found using {input_dir}')
    images.sort()
    images = np.array(images)

    output_dir = os.path.join(args.basedir, 'rerun', output_rerun, args.exptype, args.pointing[0], args.filter)
    os.makedirs(output_dir, exist_ok=True)
    logging.info(f'Writing results to {output_dir}')

    hdus = [ fits.open(image) for image in images ]

    # In debug mode just do three images or less if there aren't three
    if logging.getLogger().getEffectiveLevel() < logging.INFO:
        stride = len(hdus)/min(5,len(hdus))
        hdus = hdus[:min(3,len(hdus))]

    logging.info(f'Created HDUs for {len(hdus)} fits files from disk')

    reference_date = time.Time([mid_exposure_mjd(hdu[0]) for hdu in hdus ]).sort()
    reference_date = reference_date[len(reference_date)//2]
    logging.debug(f'Determined the reference date to be {reference_date.isot} {reference_date.mjd}')

    wcs_out, shape_out = find_optimal_celestial_wcs([hdu[1] for hdu in hdus], 
                                                    resolution=pixel_scale, projection='TAN')
    logging.info(f'Computed the optimal WCS for stacking:')
    logging.debug(f'{wcs_out}')

    if 'shiftOne' in output_dir:
        hdus = hdus[0::3]
    if 'shiftTwo' in output_dir:
        hdus = hdus[1::3]
    elif 'shiftThree' in output_dir:
        hdus = hdus[2::3]
    rates = shift_rates(1, 5, -3, 3)

    # Project the input images to the same grid using interpolation
    hdu_idx = {'data': 1, 'mask': 2, 'variance': 3}
    for rate in rates:
        stack_input = []
        logging.info(f'stacking at rate/angle set: {rate}')
        for hdu in hdus:
            data = hdu[1].data
            wcs = hdu[1].header.copy()
            dt = (mid_exposure_mjd(hdu[0]) - reference_date)
            wcs['CRVAL1'] += (rate[0]*dt*np.cos(rate[1])).to('degree').value
            wcs['CRVAL2'] += (rate[0]*dt*np.sin(rate[1])).to('degree').value 
            nddata_parts = []
            for layer in hdu_idx:
            stack_input.append((data, wcs))
        logging.info('coadd started')
        array, footprint = reproject_and_coadd(stack_input, 
                                               wcs_out, 
                                               shape_out=shape_out, 
                                               reproject_function=reproject_interp)
        logging.info('writing stack to disk')
        fits.ImageHDU(data=array, header=wcs_out.to_header()).writeto(f'{output_dir}/shifted_{ccd}_{rate[0].value}_{rate[1].value}.fits')

if __name__ == "__main__":
    main()
