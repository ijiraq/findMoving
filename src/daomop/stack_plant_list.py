from glob import glob
import logging
from matplotlib import pyplot as plt
from matplotlib import colormaps, cm
import matplotlib.transforms as transforms
from matplotlib.patches import Polygon
import numpy as np
import os
import os.path
import pickle
from scipy.optimize import curve_fit
import warnings

from astropy.coordinates import SkyCoord
from astropy import units
from astropy.io import fits
from astropy.table import Table, join, setdiff, vstack, unique, MaskedColumn
from astropy.time import Time
from astropy import units as u
from astropy import units
from astropy.wcs import WCS
from astroquery import cadc
from spherical_geometry.polygon import SphericalPolygon
from spherical_geometry.vector import lonlat_to_vector, radec_to_vector

from mp_ephem import Observation, BKOrbit

DBIMAGES = "/arc/projects/classy/dbimages"
WORKDIR = "/tmp/"
plant_list_column_names = ['id', 'ra', 'dec', 'x', 'y', 'rate', 'angle', 'dra', 'ddec', 'mag', 'amp', 'gi']


def query_cadc_for_classy_exposures(ra: float, dec: float, nights: list[int], radius: float = 0.15) -> Table:
    """
    Query the CADC for CFHT MegaCam exposures that cover the given coordinates and time.

    Args:
        ra: central right ascension, in degrees, of the classy block
        dec: central declination, in degrees, of the classy block
        nights: MJD of the UT night(s) to find exposure information for
        radius: size of the search radius around the coordinates, in degrees

    Returns:
        table of exposures that cover the given coordinates and time.
    """
    # Compute the mean of the block pointing and discovery epoch
    nights = " OR ".join([ f"INTERSECTS( INTERVAL( {Time(night).mjd}, {Time(night).mjd+1}), Plane.time_bounds ) =1 "
                           for night in nights])
    query = (
        "SELECT distinct caom2.Observation.proposal_id, "
        "caom2.Observation.target_name as target_name, "
        "min(caom2.Observation.sequenceNumber) as start_expnum, "
        "max(caom2.Observation.sequenceNumber) as end_expnum, "
        "avg(caom2.Observation.targetPosition_coordinates_cval1) as RA, "
        "avg(caom2.Observation.targetPosition_coordinates_cval2) as DEC, "
        "count(*) as nobs, "
        "SUM(caom2.Plane.time_exposure) as exptime, "
        "min(caom2.Plane.time_bounds_lower) as mjd_start, "
        "max(caom2.Plane.time_bounds_lower) as mjd_end, "
        "floor(caom2.Plane.time_bounds_upper) as night "
        "FROM caom2.Observation "
        "JOIN caom2.Plane on caom2.Observation.obsID=caom2.Plane.obsID "
        "WHERE caom2.Observation.collection='CFHT' "
        " AND caom2.Plane.time_exposure < 360 "
        " AND caom2.Plane.time_exposure > 90 "
        f" AND ({nights}) "
        f"AND INTERSECTS( CIRCLE('ICRS', {ra}, {dec}, 0.15), Plane.position_bounds ) = 1 "
        " AND caom2.Plane.calibrationLevel=1 "
        "GROUP BY caom2.Observation.proposal_id, target_name, night "
        "HAVING  SUM(caom2.Plane.time_exposure)  > 3600 "
            )

    tap_query = cadc.Cadc().create_async(query)
    tap_query.run().wait()
    tap_query.raise_if_error()
    result = tap_query.fetch_result().to_qtable()

    logging.debug(f"Query: {query} returned {len(result)} observations")

    if not len(result) > 0:
        raise ValueError(f"Query: {query} returns no observations")

    result['coord'] = SkyCoord(result['RA'], result['DEC'], unit='degree')
    result['obs_date'] = Time(result['mjd_start'], format='mjd')
    result['day_obs'] = [ x.isot[:10].replace('-','') for x in result['obs_date'] ]
    return result


def megacam_polygon(ra: float | str, dec: float | str) -> np.ndarray:
    """
    Given ra and dec return a numpy array of corners of the MegaCam field of view.
    Args:
        ra: central right ascension, in degrees or hours
        dec: central declination, in degrees

    Returns:
        polygon corners as a numpy array of shape (13, 2) with columns [ra, dec] in degrees.
    """
    coord_unit = 'degree'
    if (":" in str(ra)) | (" " in str(ra)):
        coord_unit = ('hour', 'degree')
    coordinate = SkyCoord(ra, dec, unit=coord_unit)
    corners = np.array([[+0.49220295, +0.24865226],
                        [+0.59738348, +0.24870259],
                        [+0.60126822, -0.22868328],
                        [+0.49589672, -0.23003771],
                        [+0.49237541, -0.48772592],
                        [-0.49606355, -0.49342289],
                        [-0.50266439, -0.23590507],
                        [-0.60801834, -0.23574556],
                        [-0.61018697, +0.24149839],
                        [-0.50502853, +0.24265469],
                        [-0.50048081, +0.50044066],
                        [+0.48433452, +0.50631853],
                        [+0.49220295, +0.24865226]]) * u.degree
    corners[:, 0] /= np.cos(coordinate.ra)
    corners[:, 0] += coordinate.ra
    corners[:, 1] += coordinate.dec
    return corners


def get_discovery_pointings():
    """
    Returns a dictionary of discovery pointings for the CLASSY project.
    Returns:

    """
    discovery_pointings = {
        'AS': { 'AS1': {
                'coord': SkyCoord("22:18:50","-12:15:45", unit=('hour', 'degree'))
            },
            'AS2': {
                'coord': SkyCoord("22:22:26","-11:55:56", unit=('hour', 'degree'))
            }
        },
        'MJ': {
            'MJ1': {
                'coord': SkyCoord("15:47:18","-19:07:51", unit=('hour', 'degree')),
            },
            'MJ2':  {
                'coord': SkyCoord("15:49:43","-19:17:37", unit=('hour', 'degree')),
            }
        },
        'ON': {
            'ON1': {
                'coord':  SkyCoord("03:54:58","19:36:59", unit=('hour', 'degree')),
            },
            'ON2': {
                'coord': SkyCoord("03:59:29","19:51:59", unit=('hour', 'degree')),
            }
        },
        'JF': {
            'JF1': {
                'coord': SkyCoord("08:00:24","21:17:48", unit=('hour', 'degree')),
            },
            'JF2': {
                'coord': SkyCoord("08:05:16","21:05:28", unit=('hour', 'degree')),
            }
        },
        'JA': {
            'JA1': {
                'coord': SkyCoord("20:11:28","-21:17:47", unit=('hour', 'degree')),
            },
            'JA2': {
                'coord': SkyCoord("20:12:02","-21:18:08", unit=('hour', 'degree')),
            }
        }
    }
    discovery_pointings['AS']['AS1']['nights'] = ['2022-08-22', '2022-08-23', '2022-08-26']
    discovery_pointings['AS']['AS2']['nights'] = ['2022-08-31', '2022-09-01', '2022-09-02']
    discovery_pointings['MJ']['MJ1']['nights'] = ['2023-05-16', '2023-05-21', '2023-05-23']
    discovery_pointings['MJ']['MJ2']['nights'] = ['2023-06-13', '2023-06-14', '2023-06-20']
    discovery_pointings['ON']['ON1']['nights'] = ['2022-11-24', '2022-11-25', '2022-11-26']
    discovery_pointings['ON']['ON2']['nights'] = ['2022-10-31', '2022-11-01', '2022-11-25']
    discovery_pointings['JF']['JF1']['nights'] = ['2023-01-19', '2023-01-21', '2023-01-22']
    discovery_pointings['JF']['JF2']['nights'] = ['2023-01-19', '2023-01-23', '2023-01-24']
    discovery_pointings['JA']['JA1']['nights'] = ['2023-06-16', '2023-06-19', '2023-06-22']
    discovery_pointings['JA']['JA2']['nights'] = ['2023-08-13', '2023-08-15', '2023-08-16']
    return discovery_pointings


def megacam_patch(ra: float | str, dec: float | str, edgecolor='k', lw=2, zorder=200) -> Polygon:
    """
    Given ra and dec return a pyplot Polygon that can be added to plot.

    ra and dec can be in degrees as floating points or in segisdecimal (ra in hours, dec in degrees)

    e.g.
    ra = 220.0
    dec= -15
    plt.gca().add_patch(megacam_patch(ra,dec))
    plt.xlim(ra+1,ra-1)
    plt.ylim(dec-1,dec+1)
    """
    corners = megacam_polygon(ra, dec)
    return Polygon(corners, closed=True, facecolor='None', fill='k', alpha=0.3, edgecolor=edgecolor, lw=lw,
                   zorder=zorder)


def megacam_field_contains(ra, dec, points):
    corners = megacam_polygon(ra, dec).value
    p = SphericalPolygon.from_radec(corners[:, 0], corners[:, 1], center=(ra, dec))
    inside = []
    for point in points:
        inside.append(p.contains_lonlat(point['RA'], point['DEC']))
    return np.array(inside)


def main(classy_block='AS'):
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
    logging.info("Starting CLASSY plant list stacking script")


    # Get the discovery pointings
    discovery_pointings = get_discovery_pointings()
    if classy_block not in discovery_pointings:
        raise ValueError(f"Block {classy_block} not found in discovery pointings {discovery_pointings.keys()}")

    for field in discovery_pointings[classy_block].keys():
        # Get the source coordinates and nights
        source_coord = discovery_pointings['AS'][field]['coord']
        source_nights = [ Time(night).mjd for night in discovery_pointings['AS'][field]['nights'] ]
        # Query CADC for exposures
        source_expnums = query_cadc_for_classy_exposures(source_coord.ra.deg, source_coord.dec.deg, source_nights)
    # Prepare directories
    dbimages = DBIMAGES
    source_dir = WORKDIR + 'classy_sources'
    os.makedirs(source_dir, exist_ok=True)
fake_id = None
for exposure_number in source_expnums:
    image = fits.open(f"{dbimages}/{exposure_number}/{exposure_number}p.fits")[chip_number+1]
    plant_list = Table.read(f"{dbimages}/{exposure_number}/ccd{chip_number:02d}/{exposure_number}p{chip_number:02d}.plantList",
                            format="ascii",
                            names=plant_list_column_names)
    if fake_id is None:
        row = plant_list[plant_list['mag']<22.6][3] #can change val in brackets to change source in list
        fake_id = int(row['id'])
        mpc_filename = f"{source_dir}/{object_name}/f{row['id']:.0f}.mpc"
        try:
            os.unlink(mpc_filename)
        except FileNotFoundError:
            pass
    row = plant_list[plant_list['id']==fake_id][0]
    date = Time(image.header['MJD-OBS'], format='mjd')
    with open(mpc_filename, "a") as ast:
        obs = Observation(null_observation=False,
                          note1=None,
                          note2=None,
                          minor_planet_number=f"{row['id']:0.0f}"[-5:],
                          discovery=False,
                          date=date,
                          ra=row['ra'],
                          dec=row['dec'],
                          mag=row['mag'],
                          observatory_code=568,
                          comment='injected')
        ast.write(obs.to_mpc()+'\n')