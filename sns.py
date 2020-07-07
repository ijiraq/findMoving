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

    input_dir = os.path.join(args.basedir, 'rerun', input_rerun, 'deepCoadd', args.filter, args.tract, args.patch)
    output_dir = os.path.join(args.basedir, 'rerun', output_rerun, 'deepCoadd', args.filter, args.tract, args.patch)

    images = glob.glob(input_dir + '/diff*.fits')
    images.sort()
    images = np.array(images)

    mjds = []
    for i, fn in enumerate(images):
        with fits.open(fn) as han:
            header = han[0].header

        reference_date = header['DATE-AVG']
        t = time.Time(reference_date, format='isot')
        print(reference_date)
        mjds.append(t.mjd)
    mid_mjd = np.mean(np.array(mjds))

    if 'shiftOne' in output_dir:
        images = images[0::3]
    if 'shiftTwo' in output_dir:
        images = images[1::3]
    elif 'shiftThree' in output_dir:
        images = images[2::3]

    rates = []
    for dr in np.linspace(1.0, 5.0, int((5.0 - 1.0) / 0.25) + 1):
        for dd in np.linspace(-3.0, 3.0, int(6.0 / .25) + 1):
            rates.append([dr, dd])
            logging.debug(f'Rate: {dr} and Angel:{dd}')

    """
        rates = [[4.0,-2.6],[3.5,-2.6],[3.0,-2.6],[2.5,-2.6],[2.0,-2.6],[1.5,-2.6],
        [4.0,0.0],[3.5,0.0],[3.0,0.0],[2.5,0.0],[2.0,0.0],[1.5,0.0],
             [4.0,-1.3],[3.5,-1.3],[3.0,-1.3],[2.5,-1.3],[2.0,-1.3],[1.5,-1.3],
             [4.0,1.3],[3.5,1.3],[3.0,1.3],[2.5,1.3],[2.0,1.3],[1.5,1.3]]
    """

    image_data = []
    mjds = []
    for i, fn in enumerate(images):
        with fits.open(fn) as han:
            data = han[0].data
            header = han[0].header

        if data.shape != (4100, 4100):
            logging.warning(f'Skipping {fn}! {data.shape}')
            continue

        reference_date = header['DATE-AVG']
        t = time.Time(reference_date, format='isot')
        logging.debug(f'Image {fn} taken on {reference_date}')
        mjds.append(t.mjd)
        image_data.append(np.copy(data))

    mjds = np.array(mjds)
    logging.debug(f'Loaded into memory patches taken on these dates {mjds}')
    image_data = np.array(image_data)

    # pass the list of data section into stack in order.
    args = np.argsort(mjds)
    mjds = mjds[args]
    image_data = image_data[args]

    pix_scale = abs(header['CD1_1'] * 3600.0)
    for rate in rates:
        shifted = shift(image_data, mjds, rate[0], rate[1], pix_scale,
                        rf=rep_fact, stacking_mode='MEDIAN', mid_mjd=mid_mjd)
        fits.writeto(f'{output_dir}/shifted_{rate[0]}_{rate[1]}.fits', shifted, overwrite=True)


if __name__ == "__main__":
    main()
