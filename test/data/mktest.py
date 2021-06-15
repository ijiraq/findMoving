from astropy.io import fits
from astropy.time import Time
import numpy

RERUN="rerun/diff/deepDiff/03095/HSC-R2/"

data = numpy.ones((100,100))
data[45:55,45:55]=2
data = data * 100
variance = numpy.sqrt(data)
bitmask = numpy.zeros((100,100), dtype=int)
from astropy.wcs import WCS

w = WCS({'CRPIX1': 50,'CRPIX2': 50,
         'CD1_1': 0.2/3600.0, 'CD1_2': 0, 'CD2_2': 0.2/3600.0, 'CD2_1': 0,
         'CRVAL1': 0, 'CRVAL2': 0,
         'CTYPE1': 'RA---TAN',
         'CTYPE2': 'DEC--TAN',
         'RADESYS': 'ICRS',
         'CUNIT1': 'deg',
         'CUNIT2': 'deg'})

hdulist = fits.HDUList([fits.PrimaryHDU(),
                        fits.ImageHDU(data, header=w.to_header()),
                        fits.ImageHDU(variance, header=w.to_header()),
                        fits.ImageHDU(data=bitmask, header=w.to_header())])
hdulist[0].header['EXTNAME'] = 'circumstance'
hdulist[1].header['EXTNAME'] = 'image'
hdulist[2].header['EXTNAME'] = 'variance'
hdulist[3].header['EXTNAME'] = 'mask'
                    
ccd=1
hdulist[0].header['CCD'] = f'{ccd:d}'
mjd = Time("2021-06-10T00:00:00").mjd
mjd1 = mjd
expnum = 221152
for i in range(100):
    expnum += i*2
    hdulist[0].header['EXPNUM'] = f'HSCA{expnum:07d}'
    mjd += 120/3600.0/24
    hdulist[0].header['MJD-STR']=mjd
    hdulist[0].header['MJD-END']=mjd+90/3600.0/24
    for j in range(1,3):
        hdulist[j].header['CRVAL1'] = (mjd - mjd1) * 24 * 1.0/3600.0
        hdulist[j].header['CRVAL2'] = (mjd - mjd1) * 24 * 1.0/3600.0
    hdulist.writeto(f"MASKED-{expnum:07d}-{ccd:03d}.fits", overwrite=True)
