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


def swarp(hdus, reference_idx, rate):
    """
    use the WCS to project all image to the 'reference_hdu' shifting the the CRVAL of each image by rate*dt
    :param hdus: list of HDUList
    :param reference_idx: index of reference HDUList in hdus
    :param rate: dictionary with the ra/dec shift rates.
    :return:
    """
    # Project the input images to the same grid using interpolation
    hdu_idx = {'image': 1, 'mask': 2, 'variance': 3, 'weight': 3}
    reference_date = mid_exposure_mjd(hdus[reference_idx][0])
    stack_input = []
    logging.info(f'stacking at rate/angle set: {rate}')
    ccd_data = {}
    for hdu in hdus:
        wcs_header = hdu[1].header.copy()
        dt = (mid_exposure_mjd(hdu[0]) - reference_date)
        if rate is not None:
            wcs_header['CRVAL1'] += (rate['dra'] * dt)
            wcs_header['CRVAL2'] += (rate['ddec'] * dt)
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
                                       WCS(hdus[reference_idx][1].header)))
        logging.debug(f'{stack_input[-1].header}')
    if rate is not None:
        combiner = Combiner(stack_input)
        stacked_image = combiner.average_combine()
        return fits.HDUList([fits.PrimaryHDU(header=hdus[reference_idx][0]),
                             fits.ImageHDU(data=stacked_image.data, header=hdus[reference_idx][1].header)])
    else:
        return stack_input


def shift(hdus, reference_idx, rate, rf=3, stacking_mode='MEAN'):
    """
    Orginal pixel grid expantion shift+stack code from wes.

    :rtype: fits.HDUList
    :return: combined datas after shifting at dx/dy and combined using stacking_mode.
    """
    from trippy.trippy_utils import downSample2d
    stacking_mode = STACKING_MODES.get(stacking_mode, STACKING_MODES['DEFAULT'])

    rx = rate['dra']
    ry = rate['ddec']
    logging.info(f'Shifting at ({rx},{ry})')

    # outs contains the shifted versions of the arrays after down sampling.
    outs = []
    mid_mjd = mid_exposure_mjd(hdus[reference_idx][0])
    wcs = WCS(hdus[reference_idx][1].header)
    ref_skycoord = wcs.wcs_pix2world([hdus[reference_idx][1].data.shape, ], 0)
    logging.debug(f'Reference Sky Coord {ref_skycoord}')
    logging.debug(f'Reference exposure taken at {mid_mjd.isot}')
    for hdu in hdus:
        # compute the x and y shift for image at this time and scale the size of shift for the
        # scaling factor of this shift.
        logging.debug(f'Adding exposure taken at {mid_exposure_mjd(hdu[0]).isot}')
        wcs = WCS(hdu[1].header)
        dra = (rx*(mid_exposure_mjd(hdu[0]) - mid_mjd)).decompose()
        ddec = (ry*(mid_exposure_mjd(hdu[0]) - mid_mjd)).decompose()
        logging.debug(f'Working on array {hdu[0]} of size {hdu[1].data.shape} and shifting by dx {dra} and dy {ddec}')
        # Use the WCS to determine the x/y shit to allow for different imager orientations.
        sky_coord = wcs.wcs_pix2world((hdu[1].data.shape,), 0)
        logging.debug(f'Corner of the FOV is {sky_coord}')
        # Add offset needed to aling the corner of the image with the reference image.
        dra += (ref_skycoord[0][0] - sky_coord[0][0])*units.degree
        ddec += (ref_skycoord[0][1] - sky_coord[0][1])*units.degree
        c1 = wcs.wcs_world2pix(sky_coord, 0)
        c2 = wcs.wcs_world2pix([[sky_coord[0][0]+dra.to('degree').value,
                                 sky_coord[0][1]+ddec.to('degree').value], ], 0)
        dx = int(rf*(c2[0][0]-c1[0][0]))
        dy = int(rf*(c2[0][1]-c1[0][1]))
        logging.debug(f'Translates into a up-scaled pixel shift of {dx},{dy}')
        # if hdu == hdus[reference_idx] or (dx == 0 and dy == 0):
        #    outs.append(np.array(hdu[1].data))
        #    continue
        # If this is the array that is the mid point of the stack
        # then just add to the data stack as no shifting is needed.

        # up-scale the array onto an array that is 'rf' times bigger.
        rep = np.repeat(np.repeat(hdu[1].data, rf, axis=0), rf, axis=1)
        logging.debug("Data from shape {} has been sampled into shape {}".format(hdu[1].data.shape,
                                                                                 rep.shape))
        # dest_bounds are range of the index where the data should go into
        # source_bounds are the range of the index where the data come from.
        # this creates a shift in the data, using index bounds.
        offsets = {1: dx, 0: dy}
        dest_bounds = [[0, rep.shape[0]], [0, rep.shape[1]]]
        source_bounds = [[0, rep.shape[0]], [0, rep.shape[1]]]
        for axis in offsets:
            if offsets[axis] < 0:
                dest_bounds[axis] = 0, rep.shape[axis] + offsets[axis]
                source_bounds[axis] = -offsets[axis], rep.shape[axis]
            elif offsets[axis] > 0:
                dest_bounds[axis] = offsets[axis], rep.shape[axis]
                source_bounds[axis] = 0, rep.shape[axis] - offsets[axis]
        logging.debug(f'Placing data from section {source_bounds} into {dest_bounds}')
        rep[dest_bounds[0][0]:dest_bounds[0][1], dest_bounds[1][0]:dest_bounds[1][1]] = \
            rep[source_bounds[0][0]:source_bounds[0][1], source_bounds[1][0]:source_bounds[1][1]]
        outs.append(rep)
    logging.debug(f'Stacking {len(outs)} images of shape {outs[0].shape}')
    stacked_data = stacking_mode(np.array(outs), axis=0)
    logging.debug(f'Got back stack of shape {stacked_data.shape}, downSamping...')
    stacked_data = downSample2d(stacked_data, rf)
    logging.debug(f'Down sampled image has shape {stacked_data.shape}')
    return fits.HDUList([fits.PrimaryHDU(header=hdus[reference_idx][0].header),
                         fits.ImageHDU(data=stacked_data, header=hdus[reference_idx][1].header)])


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
            dra = dr*np.cos(np.deg2rad(dd)) * units.arcsecond/units.hour
            ddec = dr*np.sin(np.deg2rad(dd)) * units.arcsecond/units.hour
            rates.append({'dra': dra, 'ddec': ddec})
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
    parser.add_argument('--swarp', action='store_true', help="Use projection to do shifts, default is pixel shifts.")
    parser.add_argument('--rectify', action='store_true', help="Rectify images to WCS of reference, otherwise "
                                                               "images must be on same grid before loading.")
    parser.add_argument('--log-level', help="What level to log at? (ERROR, INFO, DEBUG)", default="ERROR",
                        choices=['INFO', 'ERROR', 'DEBUG'])

    args = parser.parse_args()
    levels = {'INFO': logging.INFO, 'ERROR': logging.ERROR, 'DEBUG': logging.DEBUG}
    logging.basicConfig(level=levels[args.log_level])

    if args.swarp:
        stack_function = swarp
    else:
        stack_function = shift

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
    reference_idx = int(len(argsorted_hdus)//2)
    reference_filename = os.path.splitext(os.path.basename(images[reference_idx]))[0]
    logging.debug(f'Will use {reference_filename} as base name for storage.')
    logging.debug(f'Determined the reference_hdu image to be {mid_exposure_mjd(hdus[reference_idx][0]).isot}')

    if 'shiftOne' in output_dir:
        hdus = hdus[0::3]
    if 'shiftTwo' in output_dir:
        hdus = hdus[1::3]
    elif 'shiftThree' in output_dir:
        hdus = hdus[2::3]
    rates = shift_rates(1, 5, -3, 3)

    if not args.swarp and args.rectify:
        # Need to project all images to same WCS before passing to stack.
        logging.info('Swarp-ing the input images to a common projection and reference frame.')
        swarped = swarp(hdus, reference_idx, None)
        for idx in range(len(swarped)):
            hdus[idx][1].data = swarped[idx].data
            hdus[idx][1].header = swarped[idx].header
            hdus[idx][2].data = swarped[idx].mask
            hdus[idx][3].data = swarped[idx].uncertainty

    for rate in rates:
        output = stack_function(hdus, reference_idx, rate)
        logging.debug(f'Got stack result {output}')
        dx = np.round(rate['dra'].to(units.arcsecond/units.hour).value, 2)
        dy = np.round(rate['ddec'].to(units.arcsecond/units.hour).value, 2)
        output.writeto(f'{output_dir}/{reference_filename}-{dx}-{dy}.fits', overwrite=True)


if __name__ == "__main__":
    main()
