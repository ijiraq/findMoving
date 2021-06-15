"""
using lines like this 

NEW,POINTING,CCD,INDEX,X,Y,RATE,ANGLE
N,03097,005,7,761.48,1090.75,4.5,0.0

produce lines like this



#  pointing   ccd   index         x         y   rate   angle    visit                                     stack                   ra                   dec   num
#     3071     2    1297    816.73    321.14    2.5     5.0   218194   STACK-0218194-002-00-+02.50-+05.00.fits    287.1405757175446   -20.816734590737564    3

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

detections = Table.read(sys.argv[1], format='ascii')
client = vos.Client()

# visit_pattern = re.compile(r"STACK-(?P<visit>\d{6-7})-(?P<chip>\d{3})_masked-00-(?P<rate>[+-]\d{2}.\d{2})-(?P<angle>[+-]\d{2}.\d{2}).fits.fz")
visit_pattern = re.compile(r"STACK-(?P<visit>\d{6,7})-(?P<chip>\d{3})(_masked)?-00-(?P<rate>[+-]\d{2}.\d{2})-(?P<angle>[+-]\d{2}.\d{2}).fits.fz")

cutout = "[1][1:1,1:1]"
for row in detections:
    if 'NEW' in detections.colnames and row['NEW'] == 'n':
        continue
    stack_uri = f"vos:NewHorizons/S20A-OT04/STACKS_V3/{row['pointing']:05d}/{row['chip']:03d}"
    stacks = client.listdir(stack_uri)
    for stack in stacks:
        match = visit_pattern.search(stack)
        if match is None:
            sys.stderr.write(f"Failed to match {stack}\n")
            continue
        visit = match.group('visit')
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
    except:
        sys.stderr.write("FAILED for {} - {}".format(visit, row['chip']))
        sys.exit(-1)
    with fits.open(fobj) as hdu:
        ra, dec = WCS(hdu[0].header).all_pix2world(row['x'],row['y'],1)
        print("{:05d} {:03d} {:5d} {:10.2f} {:10.2f} {:10.2f} {:10.2f} {} {} {:12.6f} {:12.6f} 3".format(row['pointing'],
                                                                                                         row['chip'],
                                                                                                         row['index'],
                                                                                                         row['x'],
                                                                                                         row['y'],
                                                                                                         row['rate'],
                                                                                                         row['angle'],
                                                                                                         visit,
                                                                                                         stack_image,
                                                                                                         ra, dec))
                                                                                                
                                                                                                

