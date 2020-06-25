from astropy.io import fits
from astropy import time
import numpy as np, pylab as pyl
import glob
from trippy.trippy_utils import downSample2d


def shift(datas,mjds,rx,ry,pixScale,rf=3,stackType = 'MEDIAN', mid_mjd = None):
    """
    datas is np array of fits data
    mjds is np array of image centre times
    rx,ry are rates of motion horizontally and vertically in arcsec per hour
    pixScale is in arcsec per pix
    """

    print('Shifting at {} {} "/hr'.format(rx,ry))

    if mid_mjd is None:
        mid_mjd = np.mean(mjds)
    mid = np.argmin(np.abs(mjds-mid_mjd))

    (A,B) = datas[0].shape
    A*=rf
    B*=rf

    #dt = (mjds-mjds[mid])*24.0
    dt = (mjds - mid_mjd)*24.0
    dx = (-rx*dt/pixScale*rf+0.5).astype('int')
    dy = (-ry*dt/pixScale*rf+0.5).astype('int')

    print(dx)
    print(dy)

    outs = []
    for ii in range(len(datas)):
        print(ii+1,len(datas),dx[ii],dy[ii])
        if ii == mid:
            outs.append(datas[mid])
            continue

        #f = np.zeros((A,B),dtype = 'float64')
        rep = np.repeat(np.repeat(datas[ii],rf,axis=0),rf,axis=1)

        if dy[ii]<0:
            a,b = 0,A+dy[ii]
            aa,bb = -dy[ii],A
        elif dy[ii]>0:
            a,b = dy[ii],A
            aa,bb = 0,A-dy[ii]
        else:
            a,b = 0,A
            aa,bb = 0,A

        if dx[ii]<0:
            c,d = 0,B+dx[ii]
            cc,dd = -dx[ii],B
        elif dx[ii]>0:
            c,d = dx[ii],B
            cc,dd = 0,B-dx[ii]
        else:
            c,d = 0,B
            cc,dd = 0,B

        rep[a:b,c:d] = rep[aa:bb,cc:dd]
        outs.append(downSample2d(rep,rf))
    if stackType == 'MEDIAN':
        return np.nanmedian(outs,axis=0)
    elif stackType == 'MEAN':
        return np.nanmean(outs,axis=0)
    else:
        raise TypeException('Stack type not understood. Can be MEAN or MEDIAN only.')

if __name__ == "__main__":
    repFact = 3

    imDir = 'output'
    writeDir = 'shiftThree'
    origImDir = '/media/LSSTproc/DATA_CFHT/rerun/coadd/deepCoadd/r2/0/0,0'
    images = glob.glob(imDir+'/diff*0-0,0-24269**fits')
    images.sort()
    images = np.array(images)


    mjds = []
    for i,fn in enumerate(images):
        with fits.open(fn) as han:
            header = han[0].header

        DATE = header['DATE-AVG']
        t = time.Time(DATE,format='isot')
        print(DATE)
        mjds.append(t.mjd)
    mid_mjd = np.mean(np.array(mjds))


    if 'shiftOne' in writeDir:
        images = images[0::3]
    if 'shiftTwo' in writeDir:
        images = images[1::3]
    elif 'shiftThree' in writeDir:
        images = images[2::3]

    rates = []
    for dr in np.linspace(1.0,5.0,int((5.0-1.0)/0.25)+1):
        for dd in np.linspace(-3.0,3.0,int(6.0/.25)+1):
            rates.append([dr,dd])
            print(dr,dd)

    """
        rates = [[4.0,-2.6],[3.5,-2.6],[3.0,-2.6],[2.5,-2.6],[2.0,-2.6],[1.5,-2.6],
        [4.0,0.0],[3.5,0.0],[3.0,0.0],[2.5,0.0],[2.0,0.0],[1.5,0.0],
             [4.0,-1.3],[3.5,-1.3],[3.0,-1.3],[2.5,-1.3],[2.0,-1.3],[1.5,-1.3],
             [4.0,1.3],[3.5,1.3],[3.0,1.3],[2.5,1.3],[2.0,1.3],[1.5,1.3]]
    """

    datas = []
    mjds = []
    for i,fn in enumerate(images):
        with fits.open(fn) as han:
            data = han[0].data
            header = han[0].header

        if data.shape!=(4100,4100):
            print(f'Skipping {fn}! {data.shape}')
            continue

        #ofn = origImDir+'/'+fn.split('/diff_')[1]
        #with fits.open(ofn) as han:
        #    header = han[0].header
        #    #print(han[0].data.shape)

        DATE = header['DATE-AVG']
        t = time.Time(DATE,format='isot')
        print(DATE)
        mjds.append(t.mjd)
        #mjds.append(header['MJD-OBS'])
        #print('dont forget to add half exposure times to mjd in the line above this one.')
        #w = np.where(np.isnan(data))
        #datas[w] = -32768.0
        datas.append(np.copy(data))

    mjds = np.array(mjds)
    print(mjds)
    datas = np.array(datas)

    args = np.argsort(mjds)
    mjds = mjds[args]
    datas = datas[args]

    pixScale = abs(header['CD1_1']*3600.0)

    for rate in rates:
        shifted = shift(datas,mjds,rate[0],rate[1],pixScale,rf=repFact,stackType='MEDIAN',mid_mjd=mid_mjd)
        fits.writeto(f'{writeDir}/shifted_{rate[0]}_{rate[1]}.fits',shifted,overwrite=True)
