"""
Measure the location of a KBO and its flux in a stack of images.  Should optimize the number of images in the stack
and the rate/angle of stacking to get the best possible flux measurement. works off the 'difference' images.

We will assume the WCS of the difference images are correct.
"""
from io import BytesIO

import logging
from matplotlib import pyplot
from astropy.io import fits
from astropy.time import Time
from astropy import wcs
import numpy
from scipy.interpolate import RectBivariateSpline
import os
import vos


VOS_BASE_URI = "vos:NewHorizons/S20A-OT04"


def interped_image(data):
    """
    given an image return a function that returns the pixel values at arbitrary x/y via interpolation.

    :param data:
    :return: func
    """
    return RectBivariateSpline(range(data.shape[0]), range(data.shape[1]), data)


def shift_ra_dec(ra, dec, ra_rate, dec_rate, ref_date, this_date):
    dt = (this_date - ref_date).to('hour').value
    return ra + ra_rate*dt, dec + dec_rate*dt


def exp_date(header):
    return Time(header['DATE-AVG'], format='isot', scale=str(header['TIMESYS']).lower())


def gkern(l=5, k=5, sig=1., dx=0, dy=0):
    """\
    creates gaussian kernel with side length l and k and a sigma of sig
    """

    ax = numpy.linspace(-(l - 1) / 2., (l - 1) / 2., l)
    ay = numpy.linspace(-(k - 1) / 2., (k - 1) / 2., k)
    xx, yy = numpy.meshgrid(ax, ay)

    kernel = numpy.exp(-0.5 * (numpy.square(xx-dx) + numpy.square(yy-dy)) / numpy.square(sig))

    return kernel / numpy.sum(kernel)


def metric(images, xoffset, yoffset, xrate, yrate, ref_date):
    metric = 0
    for image in images:
        x, y = images[image]['centre']
        dx, dy = shift_ra_dec(xoffset, yoffset, xrate, yrate, ref_date, images[image]['ref_date'])
        xc = x + dx
        yc = y + dy
        k, l = images[image]['hdulist'][1].data.shape
        metric += numpy.nansum(gkern(l, k, 3.5, xc, yc) * images[image]['data'])
    return metric


def stack_inputs(images):
    d = []
    for image in images:
        d.append(images[image]['hdulist'][1].data)
    return numpy.median(d)


def get_fits_header(uri):

    filename = os.path.basename(uri)
    if os.access(filename, os.R_OK):
        header = fits.open(filename)[0].header
    else:
        c = vos.Client()
        fobj = c.open(uri, view='header')
        header = fits.Header.fromtextfile(BytesIO(fobj.read()))
        fobj.close()
    return header


def get_fits_cutout(uri, ra, dec, rad):
    filename = os.path.basename(uri)
    if os.access(filename, os.R_OK):
        hdulist = fits.open(filename)
    else:
        c = vos.Client()
        fobj = c.open(uri, view='cutout', cutout=f'CIRCLE ICRS {ra} {dec} {rad}')
        hdulist = fits.open(BytesIO(fobj.read()))
        fobj.close()
    return hdulist


def diag_plot(inputs, solution, ref_date):
    x = []
    y = []
    for filename in inputs:
        dx, dy = shift_ra_dec(solution[0], solution[1], solution[2], solution[3],
                              ref_date, inputs[filename['ref_date']])
        xc, yc = inputs[filename]['centre']
        x.append(xc + dx)
        y.append(yc + dy)
    pyplot.plot(x, y, 'o-')
    pyplot.xlabel('X-pixel')
    pyplot.xlabel('Y-pixel')
    pyplot.show()



def main():
    ra = 287.58987
    dec = -20.80109
    rad = 10/3600.0
    ra_rate = -2.41/3600.0
    dec_rate = -0.22/3600.0
    logging.basicConfig(level=logging.DEBUG,
                        format='%(relativeCreated)d %(pathname)s %(funcName)s %(lineno)d %(message)s',
                        datefmt='%I:%M:%S')
    base_dir = "HSC_May25-lsst"
    rerun = "diff"
    pointing = "03071"
    filter = "HSC-R2"
    vos_dir = os.path.join(VOS_BASE_URI, base_dir, 'rerun', rerun, pointing, filter)
    c = vos.Client()
    logging.getLogger("urllib3").setLevel(logging.DEBUG)
    filenames = c.listdir(vos_dir)
    inputs = {}
    ref_date = None
    for filename in filenames:
        logging.debug(f'Reading {filename}')
        if filename not in inputs:
            inputs[filename] = {}
        inputs[filename]['uri'] = os.path.join(vos_dir, filename)
        inputs[filename]['header'] = get_fits_header(inputs[filename]['uri'])
        inputs[filename]['ref_date'] = exp_date(inputs[filename]['header'])
        if ref_date is None:
            ref_date = inputs[filename]['ref_date']
        dt = (inputs[filename]['ref_date'] - ref_date).to('hour').value
        this_ra, this_dec = shift_ra_dec(ra, dec, ra_rate, dec_rate, ref_date, inputs[filename]['ref_date'])
        logging.debug(f'{filename} {ref_date} {dt} {ra} {dec} {this_ra} {this_dec}')
        inputs[filename]['hdulist'] = get_fits_cutout(inputs[filename]['uri'], this_ra, this_dec, rad)
        w = wcs.WCS(inputs[filename]['hdulist'][1].header)
        x, y = w.wcs_world2pix(ra, dec, 0)
        logging.debug(f'{filename} {ra},{dec}->{x},{y}')
        inputs[filename]['centre'] = x, y
        inputs[filename]['data'] = inputs[filename]['hdulist'][1].data / inputs[filename]['hdulist'][3].data
        if not os.access(filename, os.F_OK):
            inputs[filename]['hdulist'].writeto(filename)
    # fits.ImageHDU(data=stack_inputs(inputs)).writeto('stack.fits')
    peak = 0
    for xoffset in numpy.arange(-10, 11, 1):
        for yoffset in numpy.arange(-10, 11, 1):
            for xrate in numpy.arange(-10, 11, 1):
                for yrate in numpy.arange(-10, 11, 1):
                    current = metric(inputs, xoffset, yoffset, xrate, yrate, ref_date)
                    if current > peak:
                        solution = (xoffset, yoffset, xrate, yrate)
                        peak = current
                        diag_plot(inputs, solution, ref_date)
                        logging.info(f'Best fit (xo, yo, xr, yr) -> {solution}: scale: {peak}')


if __name__ == '__main__':
    main()

