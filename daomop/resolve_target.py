"""
using lines like this 

NEW,POINTING,CCD,INDEX,X,Y,RATE,ANGLE
N,03097,005,7,761.48,1090.75,4.5,0.0

produce lines like this


"""
import math
import sys
import os
from astropy.table import Table
from astropy.wcs import WCS
from astropy.io import fits
from io import BytesIO
import re
import vos
import shutil

input_filename = sys.argv[1]
output_filename = input_filename
bck_filename = f"{input_filename}.bck"
shutil.copyfile(input_filename, bck_filename)

detections = Table.read(bck_filename, format='ascii')

# Table might now have these.
if 'ra' not in detections.colnames:
    detections['ra'] = -999.9
    detections['dec'] = -999.9
    detections['mjd'] = -999.9
    detections['stack'] = -999.9

client = vos.Client()

# visit_pattern = re.compile(r"STACK-(?P<visit>\d{6-7})-(?P<chip>\d{3})_masked-00-(?P<rate>[+-]\d{2}.\d{2})-(?P<angle>[+-]\d{2}.\d{2}).fits.fz")
visit_pattern = re.compile(r"STACK-(?P<visit>\d{6,7})-(?P<chip>\d{3})(_masked)?-00-(?P<rate>[+-]\d{2}.\d{2})-(?P<angle>[+-]\d{2}.\d{2}).fits.fz")

sys.stdout.write("{:5s} {:3s} {:5s} {:10s} {:10s} {:10s} {:10s} {:7s} {:15s} {:12s} {:12s} {:2s}\n".format('pointing',
    'chip', 'id', 'x', 'y', 'rate', 'angle', 'visit', 'mjd', 'ra', 'dec', 'ntck'))
cutout = "[1][1:1,1:1]"

for row in detections:
    if 'NEW' in detections.colnames and row['NEW'] == 'n':
        continue
    # Check if previously done.
    if row['ra'] > 0:
        continue
    stack_uri = f"vos:NewHorizons/S20A-OT04/STACKS_V3/{row['pointing']:05d}/{row['chip']:03d}"
    stacks = client.listdir(stack_uri)
    for stack in stacks:
        match = visit_pattern.search(stack)
        if match is None:
            sys.stderr.write(f"Failed to match {stack}\n")
            continue
        visit = match.group('visit')
        row['visit'] = visit
        rate = float(match.group('rate'))
        angle = float(match.group('angle'))
        if math.fabs(rate - row['rate']) < 0.5 and math.fabs(angle-row['angle']) < 2.5:
            break
    if match is None:
        sys.stderr.write("FAILED to match STACK pattern")
        sys.exit(-1)
    stack_image = stack
    uri = f"{stack_uri}/{stack_image}"
    try:
        fobj = BytesIO(client.open(uri, cutout=cutout, view='cutout').read())
        with fits.open(fobj) as hdu:
            row['ra'], row['dec'] = WCS(hdu[0].header).all_pix2world(row['x'],row['y'],1)
        fobjh = BytesIO(client.open(uri, cutout="[0]", view='cutout').read())
        with fits.open(fobjh) as hhdu:
            row['mjd'] = hhdu[0].header.get('MIDMJD', hhdu[0].header.get('MJD-STR', 50000))
    except:
        sys.stderr.write("FAILED for {}".format(row))
    detections.write(output_filename, format='ascii', overwrite=True)


                                                                                                
