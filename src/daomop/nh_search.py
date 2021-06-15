import ephem
import numpy
from astropy import coordinates
from astropy import time
from astropy import units
from astropy.io import ascii
from matplotlib import pyplot
from mp_ephem import EphemerisReader as mpcread

ephem_day_0 = time.Time('1899-12-31T12:00:00', scale='utc')


def phase(r, delta, robs):
    """
    r = Sun-Object distance
    delta = Obs-Object distance
    robs = Sun-Obs distance

    taken from J-M. Petit's SurveySub.f 
    """
    denom = 2.0*r*delta
    cos = (delta**2 + r**2 - robs**2)/denom
    return numpy.arccos([x if x < 1.0 else 1.0 for x in cos]) * units.radian


def appmag(r, delta, robs, h, g=0.10):
    """
    Compute the apparent magnitude using phase and slope approach
    """
    alpha = phase(r, delta, robs)
    phi1 = numpy.exp(-3.33*(numpy.tan(alpha/2.))**0.63)
    phi2 = numpy.exp(-1.87*(numpy.tan(alpha/2.))**1.22)
    return 5*numpy.log10(r*delta/(1*units.au**2)) + h - 2.5*numpy.log10((1-g)*phi1 + g*phi2)


def build_kbos(model):
    lons = []
    lats = []
    dist = []
    kbos = []
    for m in model:
        kbo = ephem.EllipticalBody()
        kbo.name = m['name']
        kbo._a = m['a']
        kbo._e = m['e']
        kbo._epoch = 2000.0
        kbo._epoch_M = 2000.0
        kbo._M = m['M']
        kbo._om = m['omega']
        kbo._Om = m['Omega']
        kbo._inc = m['i']
        kbo._H = m['H']
        kbo._epoch = (time.Time(m['epoch'], format='mjd') - ephem_day_0).to('day').value
        kbo._epoch_M = (time.Time(m['epoch_M'], format='mjd') - ephem_day_0).to('day').value
        kbo.compute()
        lons.append(kbo.hlon)
        lats.append(kbo.hlat)
        dist.append(kbo.sun_distance)
        kbos.append(kbo)
    model['hlon'] = numpy.array(lons)*units.radian
    model['hlat'] = numpy.array(lats)*units.radian
    model['dist'] = numpy.array(dist)*units.au
    model['kbo'] = numpy.array(kbos)
    coords = coordinates.SkyCoord(model['hlon'], model['hlat'],
                                  frame='heliocentrictrueecliptic', 
                                  distance=model['dist'])
    model['coordinate'] = coords
    return


#t = ascii.read('NH_Trajectory4KBOSearch_2010-2030.csv')
t = ascii.read('sept.csv')
t['day'] = [int(x.split('/')[1]) for x in t['Time (UTC)']]
t['month'] = [int(x.split('/')[0]) for x in t['Time (UTC)']]
t['year'] = [2000 + int(x.split('/')[2]) for x in t['Time (UTC)']]
t['date_str'] = ["{:4d}-{:02d}-{:02d} 00:00:00.0".format(r['year'], r['month'], r['day']) for r in t]
t['date'] = time.Time(t['date_str'])
t['coordinate'] = coordinates.SkyCoord(x=t['x (km)'] * units.km,
                                       y=t['y (km)'] * units.km,
                                       z=t['z (km)'] * units.km,
                                       representation='cartesian', frame='icrs', obstime=t['date'])
t['n'] = 0
nh_trajectory = t

model = mpcread.load(filename='L35U.DAT')

radius = 0
plot = True

# We need to append onto ModelUsed so we get to H<13 as the limit.
# assume 'alpha = 0.4 between end of ModelUsed.dat and H=13
# so, N(H<13) = N(ModelUsed.dat) * 10**(-0.4*(min(ModelUsed.dat['H']) - 13.0))
# which is about 40 times
build_kbos(model)
new_model = model

pyplot.clf()
pyplot.hist(new_model['H'], bins=numpy.arange(5, 15, 0.1), histtype='step', cumulative=True)
pyplot.xlabel('H (mag)')
pyplot.ylabel('N (cumulative)')
pyplot.title('CKBO LF')
pyplot.yscale('log')
pyplot.savefig('input_LF.pdf', overwrite=True)
pyplot.clf()

print("# size of input model: {}".format(len(new_model)))

radius = 0*units.au
best = {}
for r in nh_trajectory:

    robs = r['Sun-NH Dist (AU)']*units.au
    if radius + 0.25*units.au > robs:
        continue
    else:
        radius = robs

    for m in new_model:
        kbo = m['kbo']
        kbo.compute((r['date']-ephem_day_0).to('day').value)
        m['hlon'] = kbo.hlon
        m['hlat'] = kbo.hlat
        m['dist'] = kbo.sun_distance
    new_model['hlon'] = numpy.array(new_model['hlon'])*units.radian
    new_model['hlat'] = numpy.array(new_model['hlat'])*units.radian
    new_model['dist'] = numpy.array(new_model['dist'])*units.au
    # sample = Table(data=(sample_idxs, hlons, hlats, dists, Hs), names=('idex', 'hlon', 'hlat', 'dist', 'H'))
    coords = coordinates.SkyCoord(new_model['hlon'],
                                  new_model['hlat'],
                                  frame='heliocentrictrueecliptic',
                                  distance=new_model['dist'])
    new_model['coordinate'] = coords
    sample = new_model

    delta = r['coordinate'].separation_3d(sample['coordinate'])
    # sample['phase'] = phase(sample['coordinate'].distance, delta, robs)

    sample['alpha'] = phase(sample['coordinate'].distance, delta, robs)
    sample['alpha'] = numpy.array(numpy.degrees(sample['alpha']))*units.degree

    nh_mag = appmag(sample['coordinate'].distance, delta, robs, sample['H'])
    gr_mag = appmag(sample['coordinate'].distance, sample['coordinate'].distance, 1.0*units.au, sample['H'])
    sample['lorri_mag'] = nh_mag
    sample['delta'] = delta
    sample['ground_mag'] = gr_mag
    sample['obs_date']=r['date'].iso

    #cond2 = numpy.all((nh_mag < 20.5, gr_mag<24.5), axis=0)
    cond2 = nh_mag < 30
    import copy
    for target in new_model[cond2]:
        if target['name'] not in best:
            best[target['name']]=copy.copy(target)
        else:
            if best[target['name']]['lorri_mag'] > target['lorri_mag']:
                best[target['name']]=copy.copy(target)

for name in best:
    target = best[name] 
    print("{}\t{:12s}\t{:5.1f}".format(target['name'], target['obs_date'], target['lorri_mag']))

