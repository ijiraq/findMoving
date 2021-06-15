import argparse
import contextlib
import logging
import os
from io import BytesIO
from .version import __version__
import numpy as np
from astropy.io import fits
from matplotlib import pyplot as pyl

# in debug mode display images using pyds9
# if pyds9 exists
try:
    from pyds9 import DS9
    DISPLAY_IMAGE=False
except ImportError:
    DISPLAY_IMAGE=False

from trippy import scamp

from . import util


def get_ref_catalog(hdulist, psf_fwhm):
    """
    Run Sextractor on a given exposure (visit/ccd) to produce a catalog of stars
    to use for masking.

    This expects that there is a file '{visit:07d}-{ccd:03d}.fits' on disk with just the
    image data (sextractor does not need to handle the MEF).

    :param hdulist: fits HDUList that holds the data to get catalog on.
    :type hdulist: fits.HDUList
    :param psf_fwhm: estimate of PSF FWHM
    :type psf_fwhm: float
    :rtype str
    :return: name of file with reference star catalog.
    """
    class MyTempFile(object):

        def __init__(self, base, fsuffix):
            self.name = f"{base}.{fsuffix}"

            if os.access(self.name, os.R_OK):
                os.unlink(self.name)

        def __del__(self):
            if logging.getLogger().getEffectiveLevel() > logging.DEBUG:
                if os.access(self.name, os.F_OK):
                    os.unlink(self.name)

    debug_mode = logging.getLogger().getEffectiveLevel() > logging.DEBUG
    _files_ = {}
    for suffix in ['cat', 'param', 'sex', 'fits', 'weight']:
        _files_[suffix] = MyTempFile("trippy_scamp", suffix)

    scamp.makeParFiles.writeConv()
    scamp.makeParFiles.writeParam(fileName=_files_['param'].name, numAps=1)
    scamp.makeParFiles.writeSex(_files_['sex'].name,
                                paramFileName=_files_['param'].name,
                                minArea=5.,
                                threshold=5.,
                                zpt=0.0,
                                aperture=psf_fwhm,
                                min_radius=2.0,
                                catalogType='FITS_LDAC',
                                saturate=60000)

    fits.writeto(filename=_files_['fits'].name, data=hdulist[1].data, header=hdulist[1].header, overwrite=True)
    fits.writeto(filename=_files_['weight'].name, data=hdulist[3].data, header=hdulist[3].header, overwrite=True)
    scamp.runSex(_files_['sex'].name, _files_['fits'].name,
                 options={'CATALOG_NAME': _files_['cat'].name,
                          'WEIGHT_IMAGE': _files_['weight'].name,
                          'WEIGHT_TYPE': 'MAP_VAR'},
                 verbose=debug_mode)
    ref_catalog = scamp.getCatalog(_files_['cat'].name,
                                   paramFile=_files_['param'].name)
    logging.debug(f"Made star catalog {ref_catalog}")

    return ref_catalog


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     fromfile_prefix_chars='@',
                                     parents=[util.base_parser])
    parser.add_argument('--psf-fwhm', help='FWHM as determined from the PSF', type=float, default=5.0)
    parser.add_argument('--clip', type=float, default=16,
                        help='Mask pixel whose variance is clip times the median variance')
    parser.add_argument('--padding-radius', help='Pad out masking by this many pixels', type=int, default=3)
    parser.add_argument('--cutout-scale', help='Box around start to examine for high variance, in psf_fwhm units.',
                        default=20, type=int)
    args = parser.parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level))

    radius_pad = args.padding_radius
    psf_fwhm = args.psf_fwhm

    show_diagnostic_plots = logging.getLogger().getEffectiveLevel() < logging.INFO

    input_rerun, output_rerun = util.parse_rerun(args.basedir, args.rerun)
    logging.debug(f"Setting input_rerun to {input_rerun}")

    corr_fn = util.get_image_list(input_rerun, exptype='CORR', visit=args.visit, ccd=args.ccd,
                                  filters=[args.filter])[0]
    diff_rerun, junk = util.parse_rerun(args.basedir, ["diff", ])
    diff_fn = util.get_image_list(diff_rerun, exptype=args.exptype, visit=args.visit,
                                  ccd=args.ccd, filters=[args.filter])[0]

    logging.debug(f'Attempting to open corr image at {corr_fn}')
    corr_hdulist = fits.open(corr_fn)
    corr = corr_hdulist[1]

    logging.debug(f'Attempting to open diff image at {diff_fn}')
    hdulist = fits.open(diff_fn)
    diff = hdulist[1]
    variance = hdulist[3]
    abs_upper_limit = args.clip * (np.nanpercentile(variance.data, 40)) ** 0.5

    if show_diagnostic_plots and DISPLAY_IMAGE:
        display = DS9('intelmask')
        display.set('zscale')
        with contextlib.closing(BytesIO()) as newFitsFile:
            fits.writeto(newFitsFile, data=diff.data, header=diff.header)
            newfits = newFitsFile.getvalue()
            display.set('fits', newfits, len(newfits))
    else:
        display = None

    
    logging.info(f'Masking {diff_fn} around bright stars found in {corr_fn} with fwhm of {psf_fwhm}')

    # cutout parameters
    cut_width = args.cutout_scale*int(psf_fwhm)
    half_cut_width = int(cut_width/2)

    # get the radial plot radii
    # Set X/Y index values at the 1/2 pixel boundaries
    # for a box that holds cut_width-X-cut_width pixels.
    x = np.arange(cut_width+1)-half_cut_width + 0.5
    y = np.arange(cut_width+1)-half_cut_width + 0.5

    inds = np.zeros((len(y), len(x), 2)).astype('int')
    for ii in range(len(y)):
        inds[ii, :, 1] = np.arange(len(x))
    for ii in range(len(x)):
        inds[:, ii, 0] = np.arange(len(y))

    coords = inds+np.array([0.5, 0.5])
    cent = np.array([cut_width/2, cut_width/2])
    r = np.sqrt(np.sum((coords-cent)**2, axis=2))
    r_reshape = r.reshape((cut_width+1)**2)

    (A, B) = corr.data.shape

    # mask_data is padded by 2*cut_width compared to the input data.
    mask_data = np.zeros((A+2*cut_width, B+2*cut_width), dtype=corr.data.dtype)
    # map the input data into the original region within mask_data.
    mask_data[cut_width:cut_width + A, cut_width:cut_width + B] = diff.data
    # mask_data is padded by 2*cut_width compared to the input data.
    vari_data = np.zeros((A+2*cut_width, B+2*cut_width), dtype=corr.data.dtype)
    # map the input data into the original region within mask_data.
    vari_data[cut_width:cut_width + A, cut_width:cut_width + B] = variance.data

    max_a = A
    max_b = B

    ref_catalog = get_ref_catalog(corr_hdulist, psf_fwhm)
    mag_bin_size = 0.5
    min_mag = max(-100, np.min(ref_catalog['MAG_AUTO']))
    max_mag = min(100, np.max(ref_catalog['MAG_AUTO']))
    mags = np.arange(max_mag, min_mag, -mag_bin_size/2)
    trim_radii = []
    trim_radius = 0


    for mag_limit in mags:

        # create a condition that pulls a magnitude range of sources from ref_catalog
        mag_select = np.where((ref_catalog['MAG_AUTO'] < mag_limit) &
                              (ref_catalog['MAG_AUTO'] > mag_limit-mag_bin_size))
        # logging.debug(f'Number of sources in mag range {mag_limit, mag_limit-mag_bin_size}: {len(mag_select[0])}')

        if len(mag_select[0]) < 10:
            trim_radii.append(trim_radius)
            continue

        vals = []
        vars = []
        for idx in mag_select[0]:
            ix, iy = int(ref_catalog['XWIN_IMAGE'][idx])+cut_width, int(ref_catalog['YWIN_IMAGE'][idx])+cut_width
            if not (0 < ix-cut_width < max_b and 0 < iy-cut_width < max_a):
                logging.warning(f'SEX gave a star outside the image {corr_fn}')
                continue
            cutout = mask_data[iy-half_cut_width:iy+half_cut_width+1, ix-half_cut_width:ix+half_cut_width+1]
            varies = vari_data[iy-half_cut_width:iy+half_cut_width+1, ix-half_cut_width:ix+half_cut_width+1]
            (a, b) = cutout.shape
            # logging.debug(f'Extracted shape of {b},{a} around {ix-half_cut_width},{iy-half_cut_width}')
            vals.append(cutout.reshape(a*b))
            vars.append(varies)
        # vals now contains a stack of x/y cutouts around the location of all ref_catalog sources in magnitude range.

        # find the pixel with a value above the abs_upper_limit cut with the largest radial distance
        # from the stellar centroid (which is the centre of the cutout box to +/- 0.5 pixels).
        try:
            med_vals = np.nanmedian(vals, axis=0)
            abs_upper_limit = args.clip*np.sqrt(np.nanmedian(vars)/len(vals))
            sky = np.nanpercentile(vals, 50)
            trim_w = np.where(np.abs(med_vals-sky) > abs_upper_limit)
            if len(trim_w[0]) > 0:
                trim_radius = np.max(r_reshape[trim_w])
            else:
                trim_radius = 0

            if show_diagnostic_plots and display is not None:
                logging.debug(f"Marking sources with {mag_limit} < MAG < {mag_limit + mag_bin_size} on diff")
                display.set('regions delete all')
                for idx in mag_select[0]:
                    display.set('regions',
                                f'image; circle({ref_catalog["XWIN_IMAGE"][idx]},'
                                f'{ref_catalog["YWIN_IMAGE"][idx]},10)')

            if show_diagnostic_plots:
                logging.debug(f"Median radial flux for {mag_limit} < MAG < {mag_limit + mag_bin_size}")
                fig = pyl.figure()
                sp = fig.add_subplot(111)
                pyl.scatter(r_reshape, med_vals-sky, alpha=0.1)
                pyl.axhline(-abs_upper_limit, color='g')
                pyl.axhline(abs_upper_limit, color='g')
                pyl.axvline(trim_radius, color='c')
                sp.set_xlabel('radial distance (pix)')
                sp.set_ylabel('pix counts (ADU)')
                pyl.ylim(-1.3*abs_upper_limit, 1.3*abs_upper_limit)
                pyl.show()
        except Exception as ex:
            logging.error(f'{corr_fn}: {ex}')
            logging.error(f'pinning the variance for {mag_limit}-{mag_limit-mag_bin_size} to 0')
            trim_radius = 0

        trim_radii.append(trim_radius)

        logging.debug(f'Trim radius for {mag_limit}:{mag_limit+mag_bin_size} -> {trim_radii[-1]}')

    num_trim_pix = 0
    for index in range(len(mags)):
        mag_range = [mags[index], mags[index]-mag_bin_size]

        if trim_radii[index] > 0:
            logging.debug(f"Trimming radius {trim_radii[index]} for mag range {mag_range}")
            trim_w = np.where(r < trim_radii[index]+radius_pad)

            w = np.where((ref_catalog['MAG_AUTO'] < mag_range[0]) & (ref_catalog['MAG_AUTO'] > mag_range[1]))
            for i in w[0]:
                ix, iy = int(ref_catalog['XWIN_IMAGE'][i])+cut_width, int(ref_catalog['YWIN_IMAGE'][i])+cut_width
                if not 0 < ix-cut_width < max_b or not 0 < iy-cut_width < max_a:
                    logging.warning(f'SEX gave a star outside the image: {corr_fn}')
                    continue
                cutout = mask_data[iy-half_cut_width:iy+half_cut_width+1, ix-half_cut_width:ix+half_cut_width+1]
                cutout[trim_w] = np.nan
                mask_data[iy-half_cut_width:iy+half_cut_width+1, ix-half_cut_width:ix+half_cut_width+1] = cutout
                num_trim_pix += len(trim_w[0])

    with fits.open(diff_fn) as hdulist:
        # pull out the part of mask_data that corresponds to the original image area.
        hdulist[1].data = mask_data[cut_width:cut_width+A, cut_width:cut_width+B]
        output_filename = diff_fn.replace(args.exptype, 'MASKED')
        hdulist[0].header['SOFTWARE'] = ( f'{__name__}-{__version__}', 'Version of daomop')
        hdulist.writeto(output_filename, overwrite=True)

    frac = num_trim_pix/(A*B)
    logging.info(f'Trimmed {num_trim_pix} or {frac} of the image area.')


if '__name__' == '__main__':
    main()
