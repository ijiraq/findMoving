import argparse
import os
import tempfile
import numpy as np
from matplotlib import pyplot as pyl
from astropy.io import fits
from trippy import scamp
import logging
from . import util


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     fromfile_prefix_chars='@',
                                     parents=[util.base_parser])
    parser.add_argument('--faint-mag-limit', help="Faintest sources to mask (in relative mags)",
                        type=float,
                        default=-8.0)
    parser.add_argument('--psf-fwhm', help='FWHM as determined from the PSF', type=float, default=5.0)
    parser.add_argument('--clip', type=int, default=16,
                        help='Mask pixel whose variance is clip times the median variance')
    parser.add_argument('--padding-radius', help='Pad out masking by this many pxiels', type=int, default=3)

    args = parser.parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level))

    faint_mag_limit = args.faint_mag_limit
    radius_pad = args.padding_radius
    psf_fwhm = args.psf_fwhm
    visit = args.visit
    ccd = args.ccd
    temp_file_list = []

    show_radial_plots = False
    if logging.getLogger().getEffectiveLevel() < logging.INFO:
        # If logging is lower than verbose (ie debug) then show plots.
        show_radial_plots = True

    input_rerun, output_rerun = util.parse_rerun(args.rerun)

    corr_dir = os.path.join(args.basedir, 'rerun', input_rerun, args.pointing, args.filter, 'corr')
    diff_dir = os.path.join(args.basedir, 'rerun', output_rerun, 'deepDiff', args.pointing, args.filter)
    logging.debug(f'Set CORR dir to: {corr_dir}')
    logging.debug(f'Set DIFF dir to: {corr_dir}')

    corr_pattern = 'CORR-{visit:07d}-{ccd:03d}.fits'
    diff_pattern = 'DIFFEXP-{visit:07d}-{ccd:03d}.fits'
    corr_fn = os.path.join(corr_dir, corr_pattern.format(visit=args.visit,
                                  ccd=args.ccd))
    diff_fn = os.path.join(diff_dir, diff_pattern.format(visit=args.visit,
                                  ccd=args.ccd))

    if not os.access(corr_fn, os.R_OK):
        # try .fz version
        corr_pattern = 'CORR-{visit:07d}-{ccd:03d}.fits.fz'
        corr_fn = os.path.join(corr_dir, corr_pattern.format(visit=args.visit,
                               ccd=args.ccd))
    logging.debug(f'Attempting to open corr image at {corr_fn}')
    with fits.open(corr_fn) as han:
        corr_data = han[1].data
        header = han[1].header

    logging.debug(f'Attempting to open diff image at {diff_fn}')
    with fits.open(diff_fn) as han:
        diff_data = han[1].data
        med_var = np.nanmedian(han[3].data)

    logging.debug(f'Determined the median variance to be: {med_var}')
    abs_upper_limit = (args.clip * med_var) ** 0.5
    logging.debug(f'Clipping abs_upper_limit is {abs_upper_limit}')
    logging.info(f'Masking {diff_fn} around bright stars found in {corr_fn} with fwhm of {psf_fwhm}')

    # cutout parameters
    cut_width = 10*int(psf_fwhm)
    half_cut_width = int(cut_width/2)

    # get the radial plot radii
    x = np.arange(cut_width+1)+0.5-5*int(psf_fwhm)
    y = np.arange(cut_width+1)+0.5-5*int(psf_fwhm)
    inds = np.zeros((len(y), len(x), 2)).astype('int')
    for ii in range(len(y)):
        inds[ii, :, 1] = np.arange(len(x))
    for ii in range(len(x)):
        inds[:, ii, 0] = np.arange(len(y))
    coords = inds+np.array([0.5, 0.5])
    cent = np.array([cut_width/2, cut_width/2])
    r = np.sqrt(np.sum((coords-cent)**2, axis=2))
    r_reshape = r.reshape((cut_width+1)**2)

    cat_filename = f'{visit:07d}-{ccd:03d}.cat'
    param_filename = f'{visit:07d}-{ccd:03d}.param'
    sex_filename = f'{visit:07d}-{ccd:03d}.sex'
    fits_filename = f'{visit:07d}-{ccd:03d}.fits'
    temp_file_list.append(cat_filename)
    temp_file_list.append(param_filename)
    temp_file_list.append(sex_filename)
    temp_file_list.append(fits_filename)
    scamp.makeParFiles.writeConv()
    scamp.makeParFiles.writeParam(fileName=param_filename, numAps=1)
    scamp.makeParFiles.writeSex(sex_filename,
                                paramFileName=param_filename,
                                minArea=5.,
                                threshold=5.,
                                zpt=0.0,
                                aperture=psf_fwhm,
                                min_radius=2.0,
                                catalogType='FITS_LDAC',
                                saturate=60000)

    fits.writeto(filename=fits_filename, data=corr_data, header=header, overwrite=True)

    scamp.runSex(sex_filename, fits_filename,
                 options={'CATALOG_NAME': cat_filename},
                 verbose=logging.getLogger().getEffectiveLevel() < logging.WARNING)
    ref_catalog = scamp.getCatalog(cat_filename,
                                   paramFile=param_filename)

    (A, B) = corr_data.shape

    exp_data = np.zeros((A+20*int(psf_fwhm), B+20*int(psf_fwhm)), dtype=corr_data.dtype)
    exp_data[10*int(psf_fwhm):10*int(psf_fwhm)+A, 10*int(psf_fwhm):10*int(psf_fwhm)+B] = diff_data
    MAX_A = A
    MAX_B = B

    mags = np.arange(faint_mag_limit, -30, -0.25)
    trim_radii = []
    trim_radius = 0
    for index in range(len(mags)):
        mag_range = [mags[index], mags[index]-0.5]
        w = np.where((ref_catalog['MAG_AUTO'] < mag_range[0]) & (ref_catalog['MAG_AUTO'] > mag_range[1]))

        if len(w[0]) < 10:
            trim_radii.append(trim_radius)
            continue

        logging.debug(f'Number of sources in mag range {mag_range}: {len(w[0])}')

        vals = []
        for i in w[0]:
            ix, iy = int(ref_catalog['XWIN_IMAGE'][i])+cut_width, int(ref_catalog['YWIN_IMAGE'][i])+cut_width
            if not 0 < ix-cut_width < MAX_B or not 0 < iy-cut_width < MAX_A+cut_width:
                logging.warning("SEX gave a star outside the image: {ref_catalog[i]}")
                continue
            cutout = exp_data[iy-half_cut_width:iy+half_cut_width+1, ix-half_cut_width:ix+half_cut_width+1]
            (a, b) = cutout.shape
            logging.debug(f'Extracting shape of {a},{b} around {ix},{iy}')
            vals.append(cutout.reshape(a*b))

        med_vals = np.nanmedian(vals, axis=0)
        trim_w = np.where(np.abs(med_vals) > abs_upper_limit)
        if len(trim_w[0]) > 0:
            trim_radius = np.max(r_reshape[trim_w])
        else:
            trim_radius = 0
        trim_radii.append(trim_radius)

        if show_radial_plots:
            fig = pyl.figure()
            sp = fig.add_subplot(111)
            pyl.scatter(r_reshape, med_vals)
            x_lim = sp.get_xlim()
            pyl.plot(x_lim, [4*med_var**0.5, abs_upper_limit])
            pyl.plot(x_lim, [-4*med_var**0.5, -abs_upper_limit])
            pyl.plot([trim_radius, trim_radius], [-500, 500])
            sp.set_ylim([-500, 500])
            sp.set_xlabel('radial distance (pix)')
            sp.set_ylabel('pix counts (ADU)')
            pyl.show()

        logging.debug(f'Trim radi list: {trim_radii}')

    num_trim_pix = 0
    for index in range(len(mags)):
        mag_range = [mags[index], mags[index]-0.5]

        if trim_radii[index] > 0:
            trim_w = np.where(r < trim_radii[index]+radius_pad)

            w = np.where((ref_catalog['MAG_AUTO'] < mag_range[0]) & (ref_catalog['MAG_AUTO'] > mag_range[1]))
            for i in w[0]:
                ix, iy = int(ref_catalog['XWIN_IMAGE'][i])+cut_width, int(ref_catalog['YWIN_IMAGE'][i])+cut_width
                if not 0 < ix-cut_width < MAX_B or not 0 < iy-cut_width < MAX_A+cut_width:
                    logging.warning("SEX gave a star outside the image: {ref_catalog[i]}")
                    continue
                cutout = exp_data[iy-half_cut_width:iy+half_cut_width+1, ix-half_cut_width:ix+half_cut_width+1]
                cutout[trim_w] = np.nan
                exp_data[iy-half_cut_width:iy+half_cut_width+1, ix-half_cut_width:ix+half_cut_width+1] = cutout
                num_trim_pix += len(trim_w[0])

    with fits.open(diff_fn) as han:
        han[1].data = exp_data[10*int(psf_fwhm):10*int(psf_fwhm)+A, 10*int(psf_fwhm):10*int(psf_fwhm)+B]
        han.writeto(diff_fn.replace('.fits', '_masked.fits'), overwrite=True)

    frac = num_trim_pix/(A*B)
    logging.info(f'Trimmed {num_trim_pix} or {frac} of the image area.')

    # Clean up temp files if this is not a DEBUG run
    if logging.getLogger().getEffectiveLevel() > logging.DEBUG:
        for filename in temp_file_list:
            if os.access(filename, os.R_OK):
                os.unlink(filename)


if '__name__' == '__main__':
    main()
