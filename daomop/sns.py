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
import sys


numpy = np

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

HSC_HDU_MAP = {'image': 1, 'mask': 2, 'variance': 3, 'weight': 3}


STACK_MASK = (2**LSST_MASK_BITS['EDGE'], 2**LSST_MASK_BITS['NO_DATA'], 2**LSST_MASK_BITS['BRIGHT_OBJECT'],
              2**LSST_MASK_BITS['SAT'], 2**LSST_MASK_BITS['INTRP'])


def mask_as_nan(data, bitmask, mask_bits=STACK_MASK):
    """
    set the mask on 'data' to include bits in mask that are set to STACK_MASK
    """
    # Check if we've been sent a mask, incase called with NDData or masked array.
    # And save that as the initial mask.
    # this is currently pretty limited.
    init_mask = None
    if isinstance(data, CCDData):
        init_mask = data.mask
        data = data.data
    # Build a boolean mask from the bits, we want masked pixels to be the one sets as ignore_flags
    # so flip_bits is True
    mask = bitfield_to_boolean_mask(bitmask, ignore_flags=mask_bits, flip_bits=True)
    if init_mask is not None:
        mask = (mask & init_mask)
    # Set masked entries to 'nan'
    data[mask] = numpy.nan
    return data


def swarp(hdus, reference_hdu, rate, hdu_idx=None, stacking_mode="MEAN"):
    """
    use the WCS to project all image to the 'reference_hdu' shifting the the CRVAL of each image by rate*dt
    :param stacking_mode: what process to use for combinnig images MEAN or MEDIAN
    :param hdu_idx: which HDU in each HDUList listed in hdus is the ImageData in?
    :param hdus: list of HDUList
    :param reference_hdu: reference HDUList in hdus
    :param rate: dictionary with the ra/dec shift rates.
    :return: fits.HDUList
    """
    # Project the input images to the same grid using interpolation
    if hdu_idx is None:
        hdu_idx = HSC_HDU_MAP
    reference_date = mid_exposure_mjd(reference_hdu[0])
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
                                       WCS(reference_hdu.header)))
        logging.debug(f'{stack_input[-1].header}')
    if rate is not None:
        combiner = Combiner(stack_input)
        if stacking_mode == 'MEDIAN':
            stacked_image = combiner.median_combine()
        else:
            stacked_image = combiner.average_combine()
        return fits.HDUList([fits.PrimaryHDU(header=reference_hdu[0]),
                             fits.ImageHDU(data=stacked_image.data, header=reference_hdu[1].header)])
    else:
        return stack_input


def shift(hdus, reference_hdu, rate, rf=3, stacking_mode='MEAN', section_size=1024):
    """
    Original pixel grid expansion shift+stack code from wes.

    :rtype: fits.HDUList
    :return: combined data after shifting at dx/dy and combined using stacking_mode.
    """
    from trippy.trippy_utils import downSample2d
    logging.info('Combining images using {stacking_mode}')
    stacking_mode = STACKING_MODES.get(stacking_mode, STACKING_MODES['DEFAULT'])

    rx = rate['dra']
    ry = rate['ddec']
    logging.info(f'Shifting at ({rx},{ry})')

    mid_mjd = mid_exposure_mjd(reference_hdu[0])
    wcs = WCS(reference_hdu[1].header)
    ref_skycoord = wcs.wcs_pix2world([reference_hdu[1].data.shape, ], 0)
    logging.debug(f'Reference Sky Coord {ref_skycoord}')
    logging.debug(f'Reference exposure taken at {mid_mjd.isot}')
    logging.info(f'Shifting {len(hdus)} to remove object motion')
    y_section_grid = np.arange(0, reference_hdu[1].data.shape[0], section_size)
    logging.debug(f'Chunk grid: y {y_section_grid}')
    x_section_grid = np.arange(0, reference_hdu[1].data.shape[1], section_size)
    logging.debug(f'Chunk grid: y {x_section_grid}')
    output_array = np.zeros(reference_hdu[1].data.shape)
    padding = 100
    for yo in y_section_grid:
        # yo,yp are the bounds were data will be inserted into output_array
        # but we need y1,y2 range of data to come from input to allow for 
        # shifting of pixel boundaries.
        yo = int(yo)
        y1 = int(max(0, yo-padding))
        yp = int(min(output_array.shape[0], yo+section_size))
        y2 = int(min(output_array.shape[0], yp+padding))
        yl = yo - y1
        yu = yl + yp - yo 
        for xo in x_section_grid:
            xo = int(xo)
            x1 = int(max(0, xo-padding))
            xp = int(min(output_array.shape[1], xo+section_size))
            x2 = int(min(output_array.shape[1], xp+padding))
            xl = xo - x1
            xu = xl + xp - xo
            logging.debug(f'Taking section {y1,y2,x1,x2} shifting, '
                          f'cutting out {yl,yu,xl,xu} '
                          f'and  placing in {yo,yp,xo,xp} ')
            # outs contains the shifted versions of the arrays after down sampling.
            outs = []
            variances = []
            for hdu in hdus:
                # compute the x and y shift for image at this time and scale the size of shift for the
                # scaling factor of this shift.
                logging.debug(f'Adding exposure taken at {mid_exposure_mjd(hdu[0]).isot}')
                wcs = WCS(hdu[1].header)
                dra = (rx*(mid_exposure_mjd(hdu[0]) - mid_mjd)).decompose()
                ddec = (ry*(mid_exposure_mjd(hdu[0]) - mid_mjd)).decompose()
                logging.debug(f'Working on array {hdu[0]} of size '
                              f'{hdu[1].data.shape} and shifting by '
                              f'dx {dra} and dy {ddec}')
                # Use the WCS to determine the x/y shit to allow for different imager orientations.
                sky_coord = wcs.wcs_pix2world((hdu[1].data.shape,), 0)
                logging.debug(f'Corner of the FOV is {sky_coord}')
                # Add offset needed to align the corner of the image with the reference image.
                dra -= (ref_skycoord[0][0] - sky_coord[0][0])*units.degree
                ddec -= (ref_skycoord[0][1] - sky_coord[0][1])*units.degree
                c1 = wcs.wcs_world2pix(sky_coord, 0)
                c2 = wcs.wcs_world2pix([[sky_coord[0][0]+dra.to('degree').value,
                                         sky_coord[0][1]+ddec.to('degree').value], ], 0)
                dx = int(rf*(c2[0][0]-c1[0][0]))
                dy = int(rf*(c2[0][1]-c1[0][1]))
                logging.debug(f'Translates into a up-scaled pixel shift of {dx},{dy}')

                # Weight by the variance. 
                rep = np.repeat(np.repeat(hdu[HSC_HDU_MAP['image']].data[y1:y2, x1:x2], rf, axis=0),
                                rf, axis=1)
                logging.debug("Data from shape {} has been sampled into shape {}".format(
                    hdu[1].data[y1:y2, x1:x2].shape, rep.shape))
                variance = np.repeat(np.repeat(hdu[HSC_HDU_MAP['variance']].data[y1:y2, x1:x2], rf, axis=0),
                                     rf, axis=1)
                rep /= variance
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
                variance[dest_bounds[0][0]:dest_bounds[0][1], dest_bounds[1][0]:dest_bounds[1][1]] = \
                    variance[source_bounds[0][0]:source_bounds[0][1], source_bounds[1][0]:source_bounds[1][1]]
                outs.append(rep)
                variances.append(variance)
            logging.debug(f'Stacking {len(outs)} images of shape {outs[0].shape}')
            logging.debug(f'Combining shifted pixels')
            stacked_variance = STACKING_MODES['MEAN'](np.array(variances), axis=0)
            stacked_data = stacking_mode(np.array(outs), overwrite_input=True, axis=0)/stacked_variance
            logging.debug(f'Got back stack of shape {stacked_data.shape}, downSampling...')
            logging.debug(f'Down sampling to original grid (poor-mans quick interp method)')
            output_array[yo:yp, xo:xp] = downSample2d(stacked_data, rf)[yl:yu, xl:xu]
    logging.debug(f'Down sampled image has shape {output_array.shape}')
    return fits.HDUList([fits.PrimaryHDU(header=reference_hdu[0].header),
                         fits.ImageHDU(data=output_array, header=reference_hdu[1].header)])


def shift_rates(r_min, r_max, r_step, angle_min, angle_max, angle_step):
    """
    @param r_min: minimum shift rate (''/hour)
    @param r_max: maximum shift rate (''/hour)
    @param angle_min: minimum angle to shift at (degrees)
    @param angle_max: maximum angle to shift at (degrees)
    @param angle_step:
    @param r_step:
    """
    rates = []
    for dd in np.linspace(angle_min, angle_max, int((angle_max - angle_min) / angle_step) + 1):
        for dr in np.linspace(r_min, r_max, int((r_max - r_min) / r_step) + 1):
            rates.append({'rate': dr, 'angle': dd})
    return rates


def mid_exposure_mjd(hdu):

    mjd_start = time.Time(hdu.header['MJD-STR'], format='mjd')
    mjd_end = time.Time(hdu.header['MJD-END'], format='mjd')
    return mjd_start + (mjd_end - mjd_start)/2.0


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('basedir', help="Root directory of LSST pipelined data")
    parser.add_argument('--pointing', help="sky patch to process (e.g. 0,0)", nargs=1)
    parser.add_argument('--rerun', help="rerun directory containing the warped difference images.", nargs=1)
    parser.add_argument('--filter', help="Filter to stack", default="HSC-R2")
    parser.add_argument('--pixel-scale', help="What should the pixel scale of the stack be? (in arc-seconds)",
                        default=0.16)
    parser.add_argument('--ccd', help="Which CCD to stack?", type=int, default=0)
    parser.add_argument('--exptype', help="What type of exposures to co-add?", default='deepDiff')
    parser.add_argument('--swarp', action='store_true', help="Use projection to do shifts, default is pixel shifts.")
    parser.add_argument('--stack-mode', choices=['MEAN', 'MEDIAN'], default='MEDIAN', help="How to combine images.")
    parser.add_argument('--rectify', action='store_true', help="Rectify images to WCS of reference, otherwise "
                                                               "images must be on same grid before loading.")
    parser.add_argument('--log-level', help="What level to log at? (ERROR, INFO, DEBUG)", default="ERROR",
                        choices=['INFO', 'ERROR', 'DEBUG'])
    parser.add_argument('--mask', action='store_true', help='set masked pixels to nan before shift/stack')
    parser.add_argument('--n-sub-stacks', default=3, help='How many sub-stacks should we produce')
    parser.add_argument('--rate-min', type=float, default=1, help='Minimum shift rate ("/hr)')
    parser.add_argument('--rate-max', type=float, default=5, help='Maximum shift rate ("/hr)')
    parser.add_argument('--rate-step', type=float, default=0.25, help='Step-size for shift rate ("/hr)')
    parser.add_argument('--angle-min', type=float, default=-3, help='Minimum angle to shift at (deg)')
    parser.add_argument('--angle-max', type=float, default=3, help='Maximum angle to shift at (deg)')
    parser.add_argument('--angle-step', type=float, default=0.25, help='Step-size for shift angle (deg)')
    parser.add_argument('--clip', type=int, default=None,
                        help='Mask pixel whose variance is clip times the median variance')
    parser.add_argument('--section-size', type=int, default=1024,
                        help='Break images into section when stacking (conserves memory)')

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
        images = images[::stride]

    hdus = [fits.open(image) for image in images]

    logging.info(f'Created HDUs for {len(hdus)} fits files from disk')

    argsorted_hdus = time.Time([mid_exposure_mjd(hdu[0]) for hdu in hdus]).argsort()
    reference_idx = int(len(argsorted_hdus)//2)
    reference_hdu = hdus[reference_idx]
    reference_filename = os.path.splitext(os.path.basename(images[reference_idx]))[0]
    logging.debug(f'Will use {reference_filename} as base name for storage.')
    logging.debug(f'Determined the reference_hdu image to be {mid_exposure_mjd(hdus[reference_idx][0]).isot}')

    # Create groups of stacks
    sub_stacks = []
    sub_images = []
    for i in range(args.n_sub_stacks):
        sub_stacks.append(hdus[i::args.n_sub_stacks])
        sub_images.append(images[i::args.n_sub_stacks])

    if not args.swarp and args.rectify:
        # Need to project all images to same WCS before passing to stack.
        logging.info('Swarp-ing the input images to a common projection and reference frame.')
        swarped = swarp(hdus, reference_hdu, None)
        for idx in range(len(swarped)):
            hdus[idx][1].data = swarped[idx].data
            hdus[idx][1].header = swarped[idx].header
            hdus[idx][2].data = swarped[idx].mask
            hdus[idx][3].data = swarped[idx].uncertainty

    if args.clip is not None:
        # Use the variance data section to mask high variance pixels from the stack.
        # mask pixels that are both high-variance AND part of a detected source.
        logging.info(f'Masking pixels in image whose variance exceeds {args.clip} times the median variance.')
        for hdu in hdus:
            hdu[HSC_HDU_MAP['variance']].header['MVAR'] = numpy.nanmedian(hdu[HSC_HDU_MAP['variance']].data)
            logging.debug(f'Median variance is {hdu[HSC_HDU_MAP["variance"]].header["MVAR"]}')
            bright_mask = hdu[HSC_HDU_MAP['variance']].data > hdu[HSC_HDU_MAP['variance']].header['MVAR']*args.clip
            detected_mask = bitfield_to_boolean_mask(hdu[HSC_HDU_MAP['mask']].data,
                                                     ignore_flags=LSST_MASK_BITS['DETECTED'],
                                                     flip_bits=True)
            logging.debug(f'Bright Mask flagged {np.sum(bright_mask)}')
            logging.debug(f'Clip setting {np.sum(bright_mask & detected_mask)} to nan')
            hdu[HSC_HDU_MAP['image']].data[bright_mask & detected_mask] = np.nan

    if args.mask:
        # set masked pixel to 'nan' before sending for stacking
        for hdu in hdus:
            hdu[HSC_HDU_MAP['image']].data = mask_as_nan(hdu[HSC_HDU_MAP['image']].data, 
                                                         hdu[HSC_HDU_MAP['mask']].data)

    for rate in shift_rates(args.rate_min, args.rate_max, args.rate_step,
                            args.angle_min, args.angle_max, args.angle_step):
        dra = rate['rate']*np.cos(np.deg2rad(rate['angle'])) * units.arcsecond/units.hour
        ddec = rate['rate']*np.sin(np.deg2rad(rate['angle'])) * units.arcsecond/units.hour
        for index, sub_stack in enumerate(sub_stacks):
            output = stack_function(sub_stack, reference_hdu, {'dra': dra, 'ddec': ddec}, 
                                    stacking_mode=args.stack_mode, section_size=args.section_size)
            logging.debug(f'Got stack result {output}')
            # Keep a history of which visits when into the stack.
            for i_index, image_name in enumerate(sub_images[index]):
                output[0].header[f'input{i_index:03d}'] = image_name
            output_filename = f'{reference_filename}-{index:02d}-{rate["rate"]:+05.2f}-{rate["angle"]:+05.2f}.fits'
            output.writeto(os.path.join(output_dir, output_filename))

    return 0


if __name__ == "__main__":
    sys.exit(main())
