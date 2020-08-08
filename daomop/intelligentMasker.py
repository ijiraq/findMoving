from astropy.io import fits
import numpy as np, pylab as pyl
import os, sys
from trippy import scamp

def main():
    faint_mag_limit = -8.0
    mask_sigma = 4.0
    radius_pad = 3
    show_radial_plots = False

    corr_fn = 'HSC_May25-lsst/rerun/processCcdOutputs/03071/HSC-R2/corr/CORR-0218190-003.fits'
    diff_fn = 'HSC_May25-lsst/rerun/diff/deepDiff/03071/HSC-R2/DIFFEXP-0218190-003.fits'
    psf_fn = 'HSC_May25-lsst/rerun/processCcdOutputs/03071/HSC-R2/corr/psfStars/CORR-0218190-003.psf.fits'

    if len(sys.argv)>1:
        corr_fn, diff_fn, psf_fn = sys.argv[1], sys.argv[2], sys.argv[3]

    with fits.open(psf_fn) as han:
        psf_fwhm = han[0].header['FWHM']

    with fits.open(corr_fn) as han:
        corr_data = han[1].data
        header = han[1].header

    with fits.open(diff_fn) as han:
        diff_data = han[1].data
        med_var = np.nanmedian(han[3].data)

    # cutout parameters
    cut_width = 10*int(psf_fwhm)
    half_cut_width = int(cut_width/2)

    # get the radial plot radii
    x = np.arange(cut_width+1)+0.5-5*int(psf_fwhm)
    y = np.arange(cut_width+1)+0.5-5*int(psf_fwhm)
    inds = np.zeros((len(y), len(x),2)).astype('int')
    for ii in range(len(y)):
        inds[ii, :, 1] = np.arange(len(x))
    for ii in range(len(x)):
        inds[:, ii, 0] = np.arange(len(y))
    coords=inds+np.array([0.5,0.5])
    cent = np.array([cut_width/2, cut_width/2])
    r=np.sqrt(np.sum((coords-cent)**2,axis=2))
    r_reshape = r.reshape((cut_width+1)**2)


    os.system('rm masker.*')
    scamp.makeParFiles.writeSex('masker.sex',
                                paramFileName='masker.param',
                                minArea=5.,
                                threshold=5.,
                                zpt=0.0,
                                aperture=psf_fwhm,
                                min_radius=2.0,
                                catalogType='FITS_LDAC',
                                saturate=60000)
    scamp.makeParFiles.writeConv()
    scamp.makeParFiles.writeParam(fileName='masker.param', numAps=1)
    junk_fn = 'junk{}.fits'.format(int((np.random.rand(1)*100000)[0]))
    fits.writeto(junk_fn, corr_data, header=header, overwrite=True)
    
    scamp.runSex('masker.sex', junk_fn ,options={'CATALOG_NAME':'masker.cat'},verbose=True)
    ref_catalog = scamp.getCatalog('masker.cat',paramFile='masker.param')
    
    os.system('rm {}'.format(junk_fn))
    
    
    (A, B) = corr_data.shape
    exp_data = np.zeros((A+20*int(psf_fwhm), B+20*int(psf_fwhm)), dtype = corr_data.dtype)
    exp_data[10*int(psf_fwhm):10*int(psf_fwhm)+A, 10*int(psf_fwhm):10*int(psf_fwhm)+B] = diff_data


    mags = np.arange(faint_mag_limit, -30, -0.25)
    trim_radii = []
    trim_radius = 0
    for l in range(len(mags)):
        mag_range = [mags[l], mags[l]-0.5]
        w = np.where((ref_catalog['MAG_AUTO']<mag_range[0]) & (ref_catalog['MAG_AUTO']>mag_range[1]))
        if len(w[0])<10:
            trim_radii.append(trim_radius)
            continue

        print(mag_range, len(w[0]))

        if show_radial_plots:
            fig = pyl.figure()
            sp = fig.add_subplot(111)

        vals = []
        for i in w[0]:
            ix, iy = int(ref_catalog['XWIN_IMAGE'][i])+cut_width, int(ref_catalog['YWIN_IMAGE'][i])+cut_width
            cutout = exp_data[iy-half_cut_width:iy+half_cut_width+1, ix-half_cut_width:ix+half_cut_width+1]

            (a, b) = cutout.shape
            vals.append(cutout.reshape(a*b))
            #pyl.scatter(r_reshape, vals[-1])

        med_vals = np.nanmedian(vals, axis=0)
        trim_w = np.where(np.abs(med_vals)>mask_sigma*med_var**0.5)
        if len(trim_w[0])>0:
            trim_radius = np.max(r_reshape[trim_w])
        else:
            trim_radius = 0
        trim_radii.append(trim_radius)

        if show_radial_plots:
            pyl.scatter(r_reshape, med_vals)
            x_lim = sp.get_xlim()
            pyl.plot(x_lim, [4*med_var**0.5, mask_sigma*med_var**0.5])
            pyl.plot(x_lim, [-4*med_var**0.5, -mask_sigma*med_var**0.5])
            pyl.plot([trim_radius, trim_radius],[-500, 500])
            sp.set_ylim([-500, 500])
            sp.set_xlabel('radial distance (pix)')
            sp.set_ylabel('pix counts (ADU)')
            pyl.show()

        print(trim_radii)

    num_trim_pix = 0
    for l in range(len(mags)):
        mag_range = [mags[l], mags[l]-0.5]

        if trim_radii[l]>0:
            trim_w = np.where(r<trim_radii[l]+radius_pad)

            w = np.where((ref_catalog['MAG_AUTO']<mag_range[0]) & (ref_catalog['MAG_AUTO']>mag_range[1]))
            for i in w[0]:
                ix, iy = int(ref_catalog['XWIN_IMAGE'][i])+cut_width, int(ref_catalog['YWIN_IMAGE'][i])+cut_width
                cutout = exp_data[iy-half_cut_width:iy+half_cut_width+1, ix-half_cut_width:ix+half_cut_width+1]
                cutout[trim_w] = np.nan
                exp_data[iy-half_cut_width:iy+half_cut_width+1, ix-half_cut_width:ix+half_cut_width+1] = cutout
                num_trim_pix += len(trim_w[0])

    with fits.open(diff_fn) as han:
        han[1].data = exp_data[10*int(psf_fwhm):10*int(psf_fwhm)+A, 10*int(psf_fwhm):10*int(psf_fwhm)+B]
        han.writeto(diff_fn.replace('.fits', '_masked.fits'), overwrite=True)

    frac = num_trim_pix/(A*B)
    print(f'Trimmed {num_trim_pix} or {frac} of the image area.')

if '__name__' == '__main__':
    main()
