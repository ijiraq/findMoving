#!/usr/bin/python3
from astropy.io import fits
import sys
import re
for filename in sys.argv[1:]:
    match = re.match(r'.*([0-9]{6})-([0-9]{3}).psf.fits', filename)
    visit=match.group(1)
    ccd=match.group(2)
    fwhm = fits.open(filename)[0].header['FWHM']
    sys.stdout.write(f'{visit} {ccd} {fwhm}\n')
