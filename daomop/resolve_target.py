"""
using lines like this 

NEW,POINTING,CCD,INDEX,X,Y,RATE,ANGLE
N,03097,005,7,761.48,1090.75,4.5,0.0

produce lines like this



#  pointing   ccd   index         x         y   rate   angle    visit                                     stack                   ra                   dec   num
#     3071     2    1297    816.73    321.14    2.5     5.0   218194   STACK-0218194-002-00-+02.50-+05.00.fits    287.1405757175446   -20.816734590737564    3

"""

import sys
import os
from astropy.table import Table
from astropy.wcs import WCS
from astropy.io import fits
from io import BytesIO
import vos

detections = Table.read(sys.argv[1], format='ascii')
visit = 221210
client = vos.Client()

cutout = "[1][1:1,1:1]"
for row in detections:
    if row['NEW'] == 'n':
        continue
    for visit in [221210 , 221208]:
        stack_image = "STACK-{:07d}-{:03d}_masked-00-{:+06.2f}-{:+06.2f}.fits.fz".format(visit, row['CCD'], 
                                                                                         row['RATE'],row['ANGLE'])
        uri = "vos:NewHorizons/S20A-OT04/STACKS_V3/{:05d}/{:03d}/{}".format(row['POINTING'], row['CCD'], stack_image)
        try:
            fobj = BytesIO(client.open(uri, cutout=cutout, view='cutout').read())
        except:
            sys.stderr.write("FAILED for {} - {}".format(visit, row['CCD']))
            continue
        break
    hdu = fits.open(fobj)
    ra, dec = WCS(hdu[0].header).all_pix2world(row['X'],row['Y'],1)
    print("{:05d} {:03d} {:5d} {:10.2f} {:10.2f} {:10.2f} {:10.2f} {} {} {:12.6f} {:12.6f} 3".format(row['POINTING'],
                                                                                                 row['CCD'],
                                                                                                 row['INDEX'],
                                                                                                 row['X'],
                                                                                                 row['Y'],
                                                                                                 -1*row['RATE'],
                                                                                                 -1*row['ANGLE'],
                                                                                                 visit,
                                                                                                 stack_image,
                                                                                                ra, dec))
                                                                                                
                                                                                                

