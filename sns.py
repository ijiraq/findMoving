import argparse
import glob
import logging
import os

import numpy as np
from astropy import time, units
from astropy.io import fits
from astropy.nddata import VarianceUncertainty, bitfield_to_boolean_mask
from astropy.wcs import WCS
from ccdproc import CCDData, wcs_project, Combiner

STACKING_MODES = {'MEDIAN': np.nanmedian,
                  'MEAN': np.nanmean,
                  'SUM': np.nansum,
                  'MAX': np.nanmax,
                  'DEFAULT': np.nanmedian}

LSST_MASK_BITS = {'BAD': 0,
                  'SAT': 1,
                  'INTRP': 2,
                  'EDGE': 4,
                  'DETECTED': 5,
                  'DETECTED_NEGATIVE': 6,
                  'SUSPECT': 7,
                  'NO_DATA': 8,
                  'CROSSTALK': 9,
                  'NOT_BLENDED': 10,
                  'UNMASKEDNAN': 11,
                  'BRIGHT_OBJECT': 12,
                  'CLIPPED': 13,
                  'INEXACT_PSF': 14,
                  'REJECTED': 15,
                  'SENSOR_EDGE': 16,
                  }

STACK_MASK = (2**LSST_MASK_BITS['EDGE'], 2**LSST_MASK_BITS['NO_DATA'], 2**LSST_MASK_BITS['BRIGHT_OBJECT'],
              2**LSST_MASK_BITS['SAT'], 2**LSST_MASK_BITS['INTRP'])


def shift(data, mjd, rx, ry, pix_scale, rf=3, stacking_mode='MEDIAN', mid_mjd=None):
    """
    Orginal pixel grid expantion shift+stack code from wes.

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
    @param r_min: minimum shift rate (''/hour)
    @param r_max: maximum shift rate (''/hour)
    @param angle_min: minimum angle to shift at (degrees)
    @param angle_max: maximum angle to shift at (degrees)
    """
    rates = []
    for dd in np.linspace(angle_min, angle_max, int((angle_max - angle_min) / .25) + 1):
        for dr in np.linspace(r_min, r_max, int((r_max - r_min) / 0.25) + 1):
            rates.append([dr * units.arcsecond / units.hour, dd * units.degree])
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

    # In debug mode just do three images or less if there aren't three
    if logging.getLogger().getEffectiveLevel() < logging.INFO:
        nimgs = min(6, len(images))
        stride = max(1, int(len(images)/nimgs-1))
        logging.debug(f'Selecting {nimgs}, every {stride} image list.')
        images = images[:nimgs:stride]

    hdus = [fits.open(image) for image in images]

    logging.info(f'Created HDUs for {len(hdus)} fits files from disk')

    argsorted_hdus = time.Time([mid_exposure_mjd(hdu[0]) for hdu in hdus]).argsort()
    reference_hdu = hdus[argsorted_hdus[len(argsorted_hdus)//2]]
    print(len(argsorted_hdus)//2)
    print(images)
    reference_filename = os.path.splitext(os.path.basename(images[argsorted_hdus[len(argsorted_hdus)//2]]))[0]
    logging.debug(f'Will use {reference_filename} as base name for storage.')
    logging.debug(f'Determined the reference_hdu image to be {mid_exposure_mjd(reference_hdu[0]).isot}')

    if 'shiftOne' in output_dir:
        hdus = hdus[0::3]
    if 'shiftTwo' in output_dir:
        hdus = hdus[1::3]
    elif 'shiftThree' in output_dir:
        hdus = hdus[2::3]
    rates = shift_rates(1, 5, -3, 3)

    # Project the input images to the same grid using interpolation
    hdu_idx = {'image': 1, 'mask': 2, 'variance': 3, 'weight': 3}
    reference_date = mid_exposure_mjd(reference_hdu[0])
    for rate in rates:
        stack_input = []
        stack_weights = []
        logging.info(f'stacking at rate/angle set: {rate}')
        ccd_data = {}
        for hdu in hdus:
            wcs_header = hdu[1].header.copy()
            dt = (mid_exposure_mjd(hdu[0]) - reference_date)
            wcs_header['CRVAL1'] += (rate[0]*dt*np.cos(rate[1].to('radian').value)).to('degree').value
            wcs_header['CRVAL2'] += (rate[0]*dt*np.sin(rate[1].to('radian').value)).to('degree').value
            for layer in hdu_idx:
                data = hdu[hdu_idx[layer]].data
                if layer == 'variance':
                    data = VarianceUncertainty(data)
                elif layer == 'mask':
                    data = bitfield_to_boolean_mask(data, ignore_flags=STACK_MASK, flip_bits=True)
                ccd_data[layer] = data
            logging.info(f'Adding {hdu[0]} to projected stack.')
            stack_input.append(wcs_project(CCDData(ccd_data['image'],
                                                   mask=ccd_data['mask'],
                                                   header=wcs_header,
                                                   wcs=WCS(wcs_header),
                                                   unit='adu',
                                                   uncertainty=ccd_data['variance']),
                                           WCS(reference_hdu[1].header)))
            # stack_weights.append(1/ccd_data['weight'].data)
        combiner = Combiner(stack_input)
        # combiner.weights = np.array(stack_weights)
        stacked_image = combiner.average_combine()
        output_hdu = fits.HDUList()
        output_hdu.append(fits.PrimaryHDU(header=reference_hdu[0].header))
        output_hdu.append(fits.ImageHDU(data=stacked_image.data, header=reference_hdu[1].header))
        output_hdu.writeto(f'{output_dir}/{reference_filename}-{rate[0].value}-{rate[1].value}.fits',
                           overwrite=True)


if __name__ == "__main__":
    main()
