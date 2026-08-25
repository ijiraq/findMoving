from unittest import TestCase
from unittest.mock import Mock
from astropy import units
from astropy.io.fits import Header
from astropy.time import Time
from . import sns
from . import util
import numpy
import mp_ephem  # noqa: F401 - registers mpc time format


class Test(TestCase):

    def test_weighted_quantile(self):
        n = numpy.random.choice(numpy.arange(50), [50, ], replace=False)
        image_stack = numpy.array([numpy.ones((10, 10))*i for i in n])
        image_weights = 0*image_stack+1/50.
        wq = sns.weighted_quantile(image_stack, 0.50001, image_weights)
        self.assertEqual(wq.shape[0], image_stack.shape[1])
        self.assertEqual(wq.shape[1], image_stack.shape[2])
        self.assertAlmostEqual(wq[5, 5], 25, 2)

    def test_position_uncertainty_pixels(self):
        orbit = Mock()
        orbit.dra = 10 * units.arcsec
        orbit.ddec = 20 * units.arcsec
        radius = sns.position_uncertainty_pixels(orbit, 0.2)
        self.assertAlmostEqual(radius, 200.0)

    def test_mid_exposure_mjd_from_date_avg(self):
        header = Header([
            ('DATE-AVG', '2022-08-20T13:29:32.725000000'),
            ('MJD-OBS', 59811.5600105),
            ('MJDEND', 59811.5635017),
        ])
        hdu = Mock(header=header)
        mid = util.mid_exposure_mjd(hdu)
        self.assertEqual(mid.scale, 'tai')
        self.assertLess(abs((mid - Time('2022-08-20T13:29:32.725000000', scale='tai')).to(units.s).value), 0.01)

    def test_mid_exposure_mjd_from_mjd_obs_end(self):
        header = Header([
            ('MJD-OBS', 59811.5600105),
            ('MJDEND', 59811.5635017),
        ])
        hdu = Mock(header=header)
        mid = util.mid_exposure_mjd(hdu)
        date_avg = Time('2022-08-20T13:29:32.725000000', scale='tai')
        self.assertLess(abs((mid - date_avg).to(units.s).value), 0.01)

    def test_mid_exposure_mpc_is_utc(self):
        header = Header([
            ('DATE-AVG', '2022-08-20T13:29:32.725000000'),
            ('MJD-OBS', 59811.5600105),
            ('MJDEND', 59811.5635017),
        ])
        hdu = Mock(header=header)
        obsdate = util.mid_exposure_mpc(hdu)
        parsed = Time(obsdate, format='mpc', scale='utc', precision=5)
        utc_mid = Time((59811.5600105 + 59811.5635017) / 2.0, format='mjd', scale='utc')
        self.assertLess(abs((parsed - utc_mid).to(units.s).value), 1.0)
