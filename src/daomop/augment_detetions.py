"""
This script takes the .detections formated file that Wes produces and augments it with RA/DEC and STACK/VISIT info.

Wes makes files with this header:

# Visit chip index    x         y    rate   angle
# pointing   ccd   index         x         y   rate   angle   visit                                             stack                   ra                   dec  nimg 


augment spits out:

# pointing ccd index x y rate angle visit stack ra dec nimg


"""
from astropy.table import Table
from astropy.wcs import WCS
from astropy.io import fits
from astropy.coordinates import SkyCoord
import sys
import argparse
import glob
from pathlib import Path

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('raw_filename', help="Filename of detection list in format Wes provides")
    parser.add_argument('aug_filename', help="Filename to store augmented version to")
    parser.add_argument('stack_dirname', help="Name of local directory containing the stacks")

    args = parser.parse_args()
    raw = Table.read(args.raw_filename, format='ascii')

    for row in raw:
        print(row)
        for path in Path(args.stack_dirname).rglob('CORR*fits'):
            w = WCS(fits.open(path)[1].header)
            ra, dec = w.all_pix2world(row['x'], row['y'], 1)
            c = SkyCoord(ra, dec, unit='degree')
            print(ra,dec,w.footprint_contains(c),path)
        sys.exit()

if __name__ == '__main__':
    main()
