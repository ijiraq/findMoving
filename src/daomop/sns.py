import argparse
import logging
import os
import tempfile
import numpy as np
import vos
from astropy import time, units
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import VarianceUncertainty, bitfield_to_boolean_mask
from astropy.wcs import WCS
from ccdproc import CCDData, wcs_project, Combiner

from . import util
from .util import get_image_list
from .version import __version__

numpy = np

previous_snapshot = None

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

HSC_HDU_MAP = {'circumstance': 0, 'image': 1, 'mask': 2, 'variance': 3, 'weight': 3}


STACK_MASK = (2**LSST_MASK_BITS['EDGE'],
              2**LSST_MASK_BITS['NO_DATA'],
              2**LSST_MASK_BITS['BRIGHT_OBJECT'],
              2**LSST_MASK_BITS['SAT'],
              2**LSST_MASK_BITS['INTRP'],
              2**LSST_MASK_BITS['REJECTED'])



def _weighted_quantile(values, quantile, sample_weight):
    """ Very close to numpy.percentile, but supports weights.
    Always overwrite=True, works on arrays with nans.

    # THIS IS NOT ACTUALLY THE PERCENTILE< BUT CLOSE ENOUGH...<

    this was taken from a stackoverflow post:
    https://stackoverflow.com/questions/21844024/weighted-percentile-using-numpy

    NOTE: quantiles should be in [0, 1]!

    :param values: numpy.array with data
    :param quantile: array-like with many quantiles needed
    :param sample_weight: array-like of the same length as `array`
    :return: numpy.array with computed quantiles.
    """
    logging.info(f'computing weighted quantile: {quantile}')
    sorter = np.argsort(values, axis=0)
    values = numpy.take_along_axis(values, sorter, axis=0)
    sample_weight = numpy.take_along_axis(sample_weight, sorter, axis=0)
    # check for inf weights, and remove
    sample_weight[numpy.isinf(sample_weight)] = 0.0
    weighted_quantiles = np.nancumsum(sample_weight, axis=0) - 0.5 * sample_weight
    weighted_quantiles /= np.nansum(sample_weight, axis=0)
    ind = np.argmin(weighted_quantiles <= quantile, axis=0)
    return np.take_along_axis(values, np.expand_dims(ind, axis=0), axis=0)[0]


STACKING_MODES['WEIGHTED_MEDIAN'] = None


def mask_as_nan(data, bitmask, mask_bits=STACK_MASK):
    """
    set the mask on 'data' to include bits in mask that are set to STACK_MASK
    """
    # Check if we've been sent a mask, in case called with NDData or masked array.
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
        mask = (mask | init_mask)
    # Set masked entries to 'nan'
    data[mask] = numpy.nan
    return data


def swarp(hdus, reference_hdu, rate, hdu_idx=None, stacking_mode="MEAN", **kwargs):
    """
    use the WCS to project all image to the 'reference_hdu' shifting the the CRVAL of each image by rate*dt
    :param stacking_mode: what process to use for combining images MEAN or MEDIAN
    :param hdu_idx: which HDU in each HDUList listed in hdus is the ImageData in?
    :param hdus: list of HDUList
    :param reference_hdu: reference HDUList in hdus
    :param rate: dictionary with the ra/dec shift rates.
    :return: fits.HDUList
    """
    # Project the input images to the same grid using interpolation
    logging.debug(f"Called with {kwargs}")
    if stacking_mode not in ['MEDIAN', 'MEAN']:
        logging.warning(f'{stacking_mode} not available for swarp stack. Setting to MEAN')
        stacking_mode = 'MEAN'
    if hdu_idx is None:
        hdu_idx = HSC_HDU_MAP
    reference_date = mid_exposure_mjd(reference_hdu[0])
    reference_header = kwargs['astheads'][reference_hdu[0].header['IMAGE']]
    reference_wcs = WCS(reference_header)
    stack_input = {}
    logging.info(f'stacking at rate/angle set: {rate}')
    ccd_data = {}
    
    for image in hdus:
        with fits.open(image) as hdu:
            wcs_header = kwargs['astheads'][hdu[0].header['IMAGE']]
            # wcs_header = hdu[1].header.copy()
            dt = (mid_exposure_mjd(hdu[0]) - reference_date)
            if rate is not None:
                wcs_header['CRVAL1'] += (rate['dra'] * dt).to('degree').value
                wcs_header['CRVAL2'] += (rate['ddec'] * dt).to('degree').value
            for layer in hdu_idx:
                data = hdu[hdu_idx[layer]].data
                if layer == 'variance':
                    data = VarianceUncertainty(data)
                elif layer == 'mask':
                    data = bitfield_to_boolean_mask(data, ignore_flags=STACK_MASK, flip_bits=True)
                ccd_data[layer] = data
            logging.info(f'Adding {hdu[0].header["IMAGE"]} to projected stack.')
            # reference_header = referece_hdu[1].header
            stack_input[image] = wcs_project(CCDData(ccd_data['image'],
                                                     mask=ccd_data['mask'],
                                                     header=wcs_header,
                                                     wcs=WCS(wcs_header),
                                                     unit='adu',
                                                     uncertainty=ccd_data['variance']),
                                             reference_wcs)

    if rate is not None:
        combiner = Combiner(stack_input.values())
        if stacking_mode == 'MEDIAN':
            stacked_image = combiner.median_combine()
        else:
            stacked_image = combiner.average_combine()
        return fits.HDUList([fits.PrimaryHDU(header=reference_hdu[0].header),
                             fits.ImageHDU(data=stacked_image.data, header=reference_header)])
    else:
        return stack_input


def down_sample_2d(inp, fr):
    """
    bin up a image by factor fr
    :param inp: input array
    :param fr: downscaling factor

    :return numpy.ndarray
    """
    new_shape = inp.shape[0]//fr, inp.shape[1]//fr
    fr = inp.shape[0]//new_shape[0], inp.shape[1]//new_shape[1]
    result = inp.reshape((new_shape[0], fr[0], new_shape[1], fr[1])).mean(axis=(-1, 1))
    return result


def up_sample_2d(output_array, input_array, rf):
    y_f = rf * input_array.shape[0]
    x_f = rf * input_array.shape[1]
    logging.debug(f'Out shape {output_array.shape} from In shape {input_array.shape} by rf')
    for y in range(0, rf):
        for x in range(0, rf):
            output_array[y:y_f:rf, x:x_f:rf] = input_array


def frameid(hdu):
    return hdu[0].header['FRAMEID']


def compute_offset(hdu, rx, ry, rf, mid_mjd, ref_skycoord):
    logging.debug(f'Compute x/y shifts from : {rx},{ry} scaled by {rf} and mid_mjd: {mid_mjd}')
    w = WCS(hdu[1].header)
    dra = (rx*(mid_exposure_mjd(hdu[0]) - mid_mjd)).decompose()
    ddec = (ry*(mid_exposure_mjd(hdu[0]) - mid_mjd)).decompose()
    logging.debug(f'Working on array {hdu[0]} of size '
                  f'{hdu[1].data.shape} and shifting by '
                  f'dx {dra} and dy {ddec}')
    sky_coord = get_centre_coord(hdu)
    # Add offset needed to align the corner of the image with the reference image.
    dra += (ref_skycoord[0] - sky_coord[0])*units.degree
    ddec += (ref_skycoord[1] - sky_coord[1])*units.degree
    try:
        _x, _y = w.wcs_world2pix(sky_coord[0], sky_coord[1], 0)
        c1 = _x, _y
        _x, _y = w.wcs_world2pix(sky_coord[0]+dra.to('degree').value,
                                 sky_coord[1]+ddec.to('degree').value, 0)
        c2 = _x, _y
        dx = (rf*(c2[0]-c1[0]))
        dy = (rf*(c2[1]-c1[1]))
    except Exception as ex:
        logging.warning(ex)
        return None, None
    logging.debug(f'Translates into a up-scaled pixel shift of {dx},{dy}')
    return dx, dy

    # Weight by the variance.
    # up_sample_2d(scaled_images[frameid(hdu)], hdu[HSC_HDU_MAP['image']].data[y1:y2, x1:x2], rf)


def get_centre_coord(hdulist):
    """
    Using the image header compute the RA/DEC of the central pixel
    :param hdulist: FITS HDUList of ImageHDU to compute central RA/DEC for.
    :type hdulist: fits.HDUList
    :return: Central RA/DEC coordinate in COOSYS of the primary WCS.
    :rtype: SkyCoord
    """
    xc = hdulist[1].data.shape[1]/2.0
    yc = hdulist[1].data.shape[0]/2.0
    sky_ra, sky_dec = WCS(hdulist[1].header).wcs_pix2world(xc, yc, 0)
    sky_coord = sky_ra, sky_dec
    logging.debug(f'Centre of the CCD is {sky_coord}')
    return sky_coord


def remap_array(image_data, dx, dy):
    """
    Offset indices of image_data by dx and dy, ie 'shift' the image array.

    :param image_data: image data that will be mapped to a new array boundary.
    :type image_data: np.array
    :type dx: float
    :type dy: float
    :return: None
    """
    # dest_bounds are range of the index where the data should go into
    # source_bounds are the range of the index where the data come from.
    # this creates a shift in the data, using index bounds.
    offsets = {1: int(dx), 0: int(dy)}
    dest_bounds = [[0, image_data.shape[0]], [0, image_data.shape[1]]]
    excluded_bounds = [[0, image_data.shape[0]], [0, image_data.shape[1]]]
    source_bounds = [[0, image_data.shape[0]], [0, image_data.shape[1]]]
    for axis in offsets:
        if offsets[axis] < 0:
            dest_bounds[axis] = 0, image_data.shape[axis] + offsets[axis]
            excluded_bounds[axis] = image_data.shape[axis] + offsets[axis], image_data.shape[axis]
            source_bounds[axis] = -offsets[axis], image_data.shape[axis]
        elif offsets[axis] >= 0:
            dest_bounds[axis] = offsets[axis], image_data.shape[axis]
            excluded_bounds[axis] = 0, offsets[axis]
            source_bounds[axis] = 0, image_data.shape[axis] - offsets[axis]
    # compute the excluded area from this image_data.
    logging.debug(f'Placing data from section {source_bounds} into {dest_bounds} and excluding {excluded_bounds}')

    image_data[dest_bounds[0][0]:dest_bounds[0][1], dest_bounds[1][0]:dest_bounds[1][1]] = \
        image_data[source_bounds[0][0]:source_bounds[0][1], source_bounds[1][0]:source_bounds[1][1]]
    image_data[excluded_bounds[0][0]:excluded_bounds[0][1], excluded_bounds[1][0]:excluded_bounds[1][1]] = np.NaN


def shift(hdus, reference_hdu, rate, rf=3, stacking_mode=None, section_size=1024, **kwargs):
    """
    Original pixel grid expansion shift+stack code from wes.

    :rtype: fits.HDUList
    :return: combined data after shifting at dx/dy and combined using stacking_mode.
    """
    logging.debug(f"Ignoring {len(kwargs)} extra arguments passed to this function.")
    if stacking_mode is None:
        stacking_mode = 'SUM'
    logging.info(f'Combining images using {stacking_mode}')
    stacking_mode = STACKING_MODES.get(stacking_mode, STACKING_MODES['DEFAULT'])

    rx = rate['dra']
    ry = rate['ddec']
    logging.debug(f'Shifting at ({rx},{ry})')

    mid_mjd = mid_exposure_mjd(reference_hdu[0])
    ref_skycoord = get_centre_coord(reference_hdu)
    logging.debug(f'Reference Sky Coord {ref_skycoord}')
    logging.debug(f'Reference exposure taken at {mid_mjd.isot}')
    logging.debug(f'Shifting {len(hdus)} to remove object motion')
    y_section_grid = np.arange(0, reference_hdu[1].data.shape[0], section_size)
    logging.debug(f'Chunk grid: y {y_section_grid}')
    x_section_grid = np.arange(0, reference_hdu[1].data.shape[1], section_size)
    logging.debug(f'Chunk grid: x {x_section_grid}')
    image_array = np.zeros(reference_hdu[1].data.shape)
    variance_array = np.zeros(reference_hdu[1].data.shape)
    # putt a padding box around our image to account for the maximum object shear
    # between the first and final image (which are about 4 hours apart)

    for yo in y_section_grid:
        # yo,yp are the bounds were data will be inserted into image_array
        # but we need y1,y2 range of data to come from input to allow for 
        # shifting of pixel boundaries.
        yo = int(yo)
        yp = int(min(image_array.shape[0], yo+section_size))
        for xo in x_section_grid:
            xo = int(xo)
            xp = int(min(image_array.shape[1], xo+section_size))
            stacks = {'image': [], 'variance': []}
            for filename in hdus:
                with fits.open(filename) as hdu_list:
                    for extension in 'image', 'variance':
                        # compute the x and y shift for image at this time and scale the size of shift for the
                        # scaling factor of this shift.
                        logging.debug(f'Adding exposure taken at {mid_exposure_mjd(hdu_list[HSC_HDU_MAP["circumstance"]]).isot} '
                                      f'reference from {mid_mjd.isot}')
                        dx, dy = compute_offset(hdu_list, rx, ry, 1, mid_mjd, ref_skycoord)
                        # x/y boundaries of cutout from this image must be shifted to account for object motion.
                        # We reduce the {x,y}1 and increase the {x,y}2 boundaries by rf to allow expanding 1/rf
                        # sub-pixel shifting.
                        logging.debug(f"Shifting by integer amount of {dx}, {dy}")
                        logging.debug(f"Lower / upper y bounds: {yo+dy-1},{yp+dy+1}")
                        logging.debug(f"Lower / upper x bounds: {xo+dx-1},{xp+dx+1}")
                        y1 = int(max(0, yo+dy-1))
                        x1 = int(max(0, xo+dx-1))
                        y2 = int(min(hdu_list[HSC_HDU_MAP[extension]].data.shape[0], yp+dy+1))
                        x2 = int(min(hdu_list[HSC_HDU_MAP[extension]].data.shape[1], xp+dx+1))
                        if y2 < y1 or x2 < x1:
                            logging.warning(f"rate shifted off Frame {filename} section {xo,yo}")
                            continue
                        logging.debug(f"Getting {y1}:{y2},{x1}:{x2} from {filename}")
                        sample_offset = ((rf * (y1 - int((yo + dy - 1))), rf * (y2 - int((yo + dy - 1)))),
                                         (rf * (x1 - int((xo + dx - 1))), rf * (x2 - int((xo + dx - 1)))))
                        logging.debug(f"Placing into {sample_offset}")
                        data = hdu_list[HSC_HDU_MAP[extension]].data[y1:y2, x1:x2]
                        logging.debug(f"Got data shape: {data.shape} from {filename}")
                        rep = np.repeat(np.repeat(data, rf, axis=0),
                                        rf, axis=1)
                        logging.debug(f"Now shape: {rep.shape} after {rf} scale up")
                        # Create a holding array of the correct size for full section at rep:
                        hold_array = np.empty((((yp-yo)+2)*rf, (xp-xo+2)*rf))
                        hold_array.fill(np.nan)
                        logging.debug(f"Created holding array of size: {hold_array.shape}")
                        hold_array[sample_offset[0][0]:sample_offset[0][1], sample_offset[1][0]:sample_offset[1][1]] = rep
                        logging.debug(f"Got {numpy.isnan(hold_array).sum()} nans and {(~numpy.isnan(hold_array)).sum()} not nan")
                        # Using the remapped array we can make the floating point part of the dx/dy shifts
                        # but still as integers, this time with rf bigger pixels.
                        dy_rep = (dy - int(dy)) * rf
                        dx_rep = (dx - int(dx)) * rf
                        logging.debug(f"tweaking by integer amount of {dx_rep},{dy_rep}")
                        # Now we have the oversampled array where does it go in our destination?
                        remap_array(hold_array, dx_rep, dy_rep)
                        logging.debug(f"remap {numpy.isnan(hold_array).sum()} nans and {(~numpy.isnan(hold_array)).sum()} not nan")
                        hold_array = hold_array[3:-3, 3:-3]
                        logging.debug(f"trimmed {numpy.isnan(hold_array).sum()} nans and {(~numpy.isnan(hold_array)).sum()} not nan")
                        stacks[extension].append(hold_array)
                        logging.debug("Data from shape {} has been sampled into shape {}".format(
                            hdu_list[HSC_HDU_MAP['image']].data[y1:y2, x1:x2].shape, hold_array.shape))

            outs = np.array(stacks['image'])
            variances = np.array(stacks['variance'])
            num_frames = numpy.sum(~numpy.isnan(outs), axis=0)
            logging.debug(f'num_frames: {num_frames[num_frames>0]}')
            logging.debug(f'Stacking {len(outs)} images of shape {outs[0].shape}')
            logging.info(f'Combining shifted pixels using {stacking_mode}')
            if stacking_mode == STACKING_MODES['WEIGHTED_MEDIAN']:
                logging.info(f'computing weighted quantile: {0.5001}')
                sorter = np.argsort(outs, axis=0)
                outs = numpy.take_along_axis(outs, sorter, axis=0)
                variances = 1./variances
                variances = numpy.take_along_axis(variances, sorter, axis=0)
                del sorter
                # check for inf weights, and remove
                variances[numpy.isinf(variances)] = 0.0
                weighted_quantiles = np.nancumsum(variances, axis=0) - 0.5 * variances
                weighted_quantiles /= np.nansum(variances, axis=0)
                ind = np.argmin(weighted_quantiles <= 0.5001, axis=0)
                del weighted_quantiles
                stacked_data = np.take_along_axis(outs, np.expand_dims(ind, axis=0), axis=0)[0]
                del outs
            else:
                stacked_data = stacking_mode(outs, axis=0)
                del outs
            print(numpy.isnan(stacked_data).sum())
            logging.info("combined finished.")
            logging.debug(f'Setting variance to mean variance / N frames')
            stacked_variance = STACKING_MODES['MEAN'](variances, axis=0)/num_frames
            print(numpy.isnan(stacked_variance).sum())
            del variances
            logging.debug(f'Got back stack of shape {stacked_data.shape}, downSampling...')
            logging.debug(f'Down sampling to original grid (poor-mans quick interp method)')
            image_array[yo:yp, xo:xp] = down_sample_2d(stacked_data, rf)
            variance_array[yo:yp, xo:xp] = down_sample_2d(stacked_variance, rf)
    logging.debug(f'Down sampled image has shape {image_array.shape}')
    hdu_list = fits.HDUList([fits.PrimaryHDU(header=reference_hdu[0].header),
                             fits.ImageHDU(data=image_array, header=reference_hdu[HSC_HDU_MAP['image']].header),
                             fits.ImageHDU(data=variance_array, header=reference_hdu[HSC_HDU_MAP['variance']].header)])
    hdu_list[1].header['EXTNAME'] = 'STACK'
    hdu_list[2].header['EXTNAME'] = 'VARIANCE'
    return hdu_list


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
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     fromfile_prefix_chars='@',
                                     parents=[util.base_parser])

    parser.add_argument('--pixel-scale', help="What should the pixel scale of the stack be? (in arc-seconds)",
                        default=0.16)
    parser.add_argument('--swarp', action='store_true', help="Use projection to do shifts, default is pixel shifts.")
    parser.add_argument('--stack-mode', choices=STACKING_MODES.keys(),
                        default='WEIGHTED_MEDIAN', help="How to combine images.")
    parser.add_argument('--rectify', action='store_true', help="Rectify images to WCS of reference, otherwise "
                                                               "images must be on same grid before loading.")
    parser.add_argument('--mask', action='store_true', help='set masked pixels to nan before shift/stack')
    parser.add_argument('--n-sub-stacks', default=3, type=int, help='How many sub-stacks should we produce')
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
    parser.add_argument('--centre', default=None, help="only stack data around this RA/DEC (decimal degree) centre",
                        nargs=2, type=float)
    parser.add_argument('--group', action='store_true', help='Make stacks time grouped instead of striding.')

    args = parser.parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level))

    reruns = args.rerun[0].split(":")
    if len(reruns) > 2:
        raise ValueError("Don't know what to do with more then 2 rerun directories.")

    input_rerun, output_rerun = util.parse_rerun(args.basedir, args.rerun)

    output_dir = os.path.join(output_rerun, args.exptype, args.pointing, args.filter)
    os.makedirs(output_dir, exist_ok=True)
    logging.info(f'Writing results to {output_dir}')

    if args.swarp:
        stack_function = swarp
    else:
        stack_function = shift

    rates = shift_rates(args.rate_min, args.rate_max, args.rate_step,
                        args.angle_min, args.angle_max, args.angle_step)
    logging.info(f'Shift-and-Stacking the following list of rate/angle pairs: '
                 f'{[(rate["rate"],rate["angle"]) for rate in rates]}')

    logging.info(f'Loading all images matching pattern: {input_rerun}')
    images = np.array(get_image_list(input_rerun, args.exptype, ccd=args.ccd,
                                     visit=args.visit, filters=[args.pointing, args.filter]))

    # check if there are astrometric headers to overwrite the WCS in the HDU
    ast_path = os.path.join(input_rerun.replace('processCcdOutputs', 'diff').replace('diff', 'asthead'), 'mega')
    astheads = {}
    for image in images:
        filename = os.path.basename(image)
        ast_filename = os.path.join(ast_path, 
                                    filename.replace('DIFFEXP', 'CORR').replace('.fits', '.mega.head'))
        if os.access(ast_filename, os.R_OK):
            astheads[filename] = fits.Header.fromtextfile(ast_filename)
        else:
            with fits.open(image) as hdu_list:
                astheads[filename] = hdu_list[1].header
            logging.warning(f"Failed to get astheader for {ast_filename} using {filename}")
        
    if not len(images) > 0:
        raise OSError(f'No images found using {input_rerun}')

    # Organize images in MJD order.
    mjds = []
    logging.info(f"Sorting list of {len(images)} based on mjd")
    for image in images:
        try:
            with fits.open(image) as hdu:
                mjds.append(time.Time(mid_exposure_mjd((hdu[0]))))
        except Exception as ex:
            logging.error(str(ex))
            logging.error(f"Failed to open {image}")
            del images[image]

    ind = np.argsort(mjds)
    # sort the images by mjd
    images = images[ind]

    # In debug mode just do three images or less if there aren't three
    if logging.getLogger().getEffectiveLevel() < logging.INFO:
        num_of_images = min(6, len(images))
        stride = max(1, int(len(images)/num_of_images-1))
        logging.debug(f'Selecting every {stride}th images, for total of {num_of_images}')
        images = images[::stride]

    # do the stacking in groups of images as set from the CL.
    reference_hdu = None
    for index in range(args.n_sub_stacks):
        if not args.group:
            # stride the image list
            sub_images = images[index::args.n_sub_stacks]
            reference_idx = int(len(images) // 2)
            reference_image = images[reference_idx]
            reference_hdu = fits.open(images[reference_idx])
            reference_hdu[0].header['IMAGE'] = os.path.basename(reference_image)
            reference_filename = os.path.splitext(os.path.basename(images[reference_idx]))[0][8:]
        else:
            # group images by time
            start_idx = len(images)//args.n_sub_stacks*index
            start_idx = int(max(0, start_idx))
            end_idx = len(images)//args.n_sub_stacks*(index+1)
            end_idx = int(min(len(images), end_idx))
            sub_images = images[start_idx:end_idx]

        hdus = {}
        for image in sub_images:
            with fits.open(image, mode='update') as hdulist:
                hdulist[0].header['IMAGE'] = os.path.basename(image)
            hdus[image] = os.path.basename(image)

        if args.group:
            reference_idx = int(len(sub_images) // 2)
            reference_hdu = fits.open(sub_images[reference_idx])

        if not args.swarp and args.rectify:
            # Need to project all images to same WCS before passing to stack.
            logging.info('Swarp-ing the input images to a common projection and reference frame.')
            swarped_hdus = swarp(hdus, reference_hdu, None)
            for idx in swarped_hdus:
                # This HDUList order of 1 for science data and 2 for mask and 3 for uncertainty
                # is a hard coding of the LSST Pipeline order from 2020.
                image = swarped_hdus[idx]
                hdus[idx][1].header['RECTIFY'] = 'DONE'
                hdus[idx][1].data = image.data
                hdus[idx][1].header = image.header
                hdus[idx][2].data = image.mask
                hdus[idx][3].data = image.uncertainty
                hdus.close()

        if args.centre is not None:
            centre = SkyCoord(args.centre[0], args.centre[1], unit='degree')
            box_size = args.section_size//2
            logging.info(f'Extracting box of 1/2 width {box_size} pixels around {centre}')
            chdus = {}
            for key in hdus:
                with fits.open(key) as hdul:
                    image = hdul[0].header['IMAGE']
                    w = WCS(astheads[image])
                    try:
                        x, y = w.all_world2pix(args.centre[0], args.centre[1], 0)
                    except Exception as ex:
                        logging.error(str(ex))
                        hdul = None
                        continue
                    logging.info(f'{image} {centre.to_string(style="hmsdms", sep=":")} -> {x},{y}')
                    x1 = int(max(0, x - box_size))
                    x2 = int(min(hdul[1].header['NAXIS1'], x1 + 2*box_size))
                    x1 = int(max(0, x2 - 2*box_size))
                    y1 = int(max(0, y - box_size))
                    y2 = int(min(hdul[1].header['NAXIS2'], y1 + 2*box_size))
                    y1 = int(max(0, y2 - 2*box_size))
                    logging.info(f'{hdul[0].header["FRAMEID"]} -> [{y1}:{y2},{x1}:{x2}]')
                    for idx in range(1, 4):
                        try:
                            data = hdul[idx].data[y1:y2, x1:x2]
                        except Exception as ex:
                            logging.error(str(ex))
                            logging.error(f"Extracting [{y1}:{y2},{x1}:{x2}] from {hdul[0].header['FRAMEID']}")
                            data = np.ones((y2-y1+1)*(x2-x1+1))*np.nan
                            data.shape = y2-y1+1, x2-x1+1
                        hdul[idx].data = data
                        hdul[idx].header['XOFFSET'] = x1
                        hdul[idx].header['YOFFSET'] = y1
                    astheads[image]['CRPIX1'] -= x1
                    astheads[image]['CRPIX2'] -= y1
                    tt_file = tempfile.NamedTemporaryFile()
                    hdul.writeto(tt_file)
                    chdus[tt_file.name] = hdul

                    hdul.close()
            hdus = chdus

        hdus2 = {}
        # Remove from the list HDULists where there is no data left (after cutout)
        for key in hdus:
            if hdus[key] is not None:
                hdus2[key] = hdus[key]
        hdus = hdus2

        if args.clip is not None:
            # Use the variance data section to mask high variance pixels from the stack.
            # mask pixels that are both high-variance AND part of a detected source.
            logging.info(f'Masking pixels in image whose variance exceeds {args.clip} times the median variance.')
            for key in hdus:
                with fits.open(key, mode='update') as hdu:
                    if hdu[HSC_HDU_MAP['image']].header.get('CLIP', None) is not None:
                        continue
                    hdu[HSC_HDU_MAP['image']].header['CLIP'] = args.clip
                    hdu.flush()
                    hdu[HSC_HDU_MAP['variance']].header['MVAR'] = (numpy.nanmedian(hdu[HSC_HDU_MAP['variance']].data),
                                                                   'Median variance')
                    logging.debug(f'Median variance is {hdu[HSC_HDU_MAP["variance"]].header["MVAR"]}')
                    bright_mask = hdu[HSC_HDU_MAP['variance']].data > hdu[HSC_HDU_MAP['variance']].header['MVAR']*args.clip
                    detected_mask = bitfield_to_boolean_mask(hdu[HSC_HDU_MAP['mask']].data,
                                                             ignore_flags=LSST_MASK_BITS['DETECTED'],
                                                             flip_bits=True)
                    logging.debug(f'Bright Mask flagged {np.sum(bright_mask)}')
                    hdu[HSC_HDU_MAP['image']].data[bright_mask & detected_mask] = np.nan
                    logging.debug(f'Clip setting {np.sum(bright_mask & detected_mask)} to nan')
                    hdu[HSC_HDU_MAP['variance']].data[bright_mask & detected_mask] = np.nan
                    # update the local version with the clipped data, but not if this is a cutout.

        if args.mask:
            # set masked pixel to 'nan' before sending for stacking
            for key in hdus:
                with fits.open(key, mode='update') as hdu:
                    if not hdu[HSC_HDU_MAP['image']].header.get('MASKED', False):
                        hdu[HSC_HDU_MAP['image']].data = mask_as_nan(hdu[HSC_HDU_MAP['image']].data,
                                                                     hdu[HSC_HDU_MAP['mask']].data)
                        hdu[HSC_HDU_MAP['variance']].data = mask_as_nan(hdu[HSC_HDU_MAP['variance']].data,
                                                                        hdu[HSC_HDU_MAP['mask']].data)
                        hdu[HSC_HDU_MAP['image']].header['MASKED'] = True

        for rate in rates:
            dra = rate['rate']*np.cos(np.deg2rad(rate['angle'])) * units.arcsecond/units.hour
            ddec = rate['rate']*np.sin(np.deg2rad(rate['angle'])) * units.arcsecond/units.hour
            if args.group:
                int_rate = int(rate["rate"]*10)
                int_angle = int((rate['angle'] % 360)*10)
                expnum = f'{int(args.pointing)}{int_rate:02d}{int_angle:04d}{index}'
                output_filename = f'{expnum}p{args.ccd:02d}.fits'
            else:
                expnum = reference_hdu[0].header.get('FRAMEID', 'HSCA0000000').replace('HSCA', '')
                output_filename = f'STACK-{reference_filename}-{index:02d}-' \
                                  f'{rate["rate"]:+06.2f}-{rate["angle"]:+06.2f}.fits.fz'
            output_uri = f"vos:NewHorizons/S20A-OT04/STACKS_V4/{int(args.pointing):05d}/{int(args.ccd):03d}"
            if vos.Client().access(f"{output_uri}/{output_filename}"):
                logging.warning(f"{output_filename} already in {output_uri}, skipping")
                continue
            output_filename = os.path.join(output_dir, output_filename)
            if os.access(output_filename, os.R_OK):
                logging.warning(f'{output_filename} exists, skipping')
                continue
            output = stack_function(hdus, reference_hdu, {'dra': dra, 'ddec': ddec},
                                    stacking_mode=args.stack_mode, section_size=args.section_size, astheads=astheads)
            logging.debug(f'Got stack result {output}, writing to {output_filename}')
            # Keep a history of which visits when into the stack.
            output[0].header['SOFTWARE'] = f'{__name__}-{__version__}'
            output[0].header['NCOMBINE'] = (len(hdus), 'Number combined')
            output[0].header['COMBALGO'] = (args.stack_mode, 'Stacking mode')
            output[0].header['RATE'] = (rate['rate'], 'arc-second/hour')
            output[0].header['ANGLE'] = (rate['angle'], 'degree')
            output[0].header['DRA'] = (dra.value, str(dra.unit))
            output[0].header['DDEC'] = (ddec.value, str(ddec.unit))
            output[0].header['CCDNUM'] = (args.ccd, 'CCD NUMBER or DETSER')
            output[0].header['EXPNUM'] = (expnum, '[int(pointing)][rate*10][(angle%360)*10][index]')
            output[0].header['MIDMJD'] = (mid_exposure_mjd(output[0]).mjd, " MID Exposure")
            output[0].header['ASTLEVEL'] = 1
            for i_index, image_name in enumerate(sub_images):
                output[0].header[f'input{i_index:03d}'] = os.path.basename(image_name)
            output.writeto(output_filename)

    return 0


if __name__ == "__main__":
    import sys
    status = main()
    sys.exit(status)
