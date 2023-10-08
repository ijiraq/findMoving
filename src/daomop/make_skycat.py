from astropy.io import fits
from astropy.time import Time
from astropy.table import Table
import sys
from astropy.wcs import WCS
import numpy, logging
import re

rows = []
for filename in sys.argv[1:]:
    print(filename)
    match = re.search('(?P<pointing>[0-9]{4})/.*/CORR-(?P<visit>[0-9]{7})-(?P<ccd>[0-9]{3})', filename)
    if  match is None:
        logging.info(f"Skipping {filename}")
        continue
    
    with fits.open(filename) as hdulist:
        footprint = numpy.array(WCS(hdulist[1].header).calc_footprint()).T
        mjdobs = Time(hdulist[0].header['DATE-AVG'], scale='tai').mjd
    ramin = footprint[0].min()
    ramax = footprint[0].max()
    demin = footprint[1].min()
    demax = footprint[1].max()
    
    rows.append([match.group('pointing'), mjdobs, match.group('ccd'), ramin, ramax, demin, demax])

rows = numpy.array(rows)
# pointing mjdobs ccd ramin ramax decmin decmax
Table(rows, names=['pointing', 'mjdobs', 'ccd', 'ramin', 'ramax', 'decmin', 'decmax']).write(f'{match.group("visit")}_skycat.txt', format='ascii.commented_header')
