from unittest import TestCase
from unittest.mock import Mock
from astropy import units
from astropy.coordinates import SkyCoord
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

    def test_orbit_rate_grid_odd(self):
        rates = sns.orbit_rate_grid(2.0, 45.0, 5, 0.1)
        self.assertEqual(len(rates), 5)
        self.assertEqual([r['rate'] for r in rates], [1.8, 1.9, 2.0, 2.1, 2.2])
        self.assertEqual(rates[2]['rate'], 2.0)
        self.assertTrue(all(r['angle'] == 45.0 for r in rates))

    def test_orbit_rate_grid_even(self):
        rates = sns.orbit_rate_grid(2.0, 10.0, 4, 0.1)
        self.assertEqual(len(rates), 4)
        self.assertEqual([r['rate'] for r in rates], [1.8, 1.9, 2.0, 2.1])
        self.assertEqual(rates[2]['rate'], 2.0)

    def test_get_centre_uses_reference_epoch(self):
        header = Header([('DATE-AVG', '2022-08-20T13:29:32.725000000')])
        hdu0 = Mock(header=header)
        reference = [hdu0]
        expected = SkyCoord(150.0, -10.0, unit='degree')
        orbit = Mock()
        orbit.coordinate = expected
        centre = sns.get_centre(reference, orbit)
        orbit.predict.assert_called_once()
        called_time = orbit.predict.call_args[0][0]
        self.assertEqual(called_time.isot, '2022-08-20T13:29:32.725')
        self.assertIs(centre, expected)

    def test_exposure_id_from_diffexp_filename(self):
        self.assertEqual(util.exposure_id_from_name('DIFFEXP-2779494-00.fits'), 2779494)
        self.assertEqual(util.exposure_id_from_name('DIFFEXP-2779494-2779500-00.fits'), 2779494)
        self.assertEqual(util.exposure_id_from_name('/path/calexp-02779494-012.fits'), 2779494)

    def test_exposure_id_from_filename_keyword(self):
        self.assertEqual(util.exposure_id_from_filename_keyword('2779465o'), 2779465)
        self.assertEqual(util.exposure_id_from_filename_keyword('2779494O'), 2779494)
        self.assertEqual(util.exposure_id_from_filename_keyword('2779494'), 2779494)

    def test_frame_label_prefers_expnum_then_filename(self):
        header = Header([('EXPNUM', 2779494), ('FILENAME', '9999999o')])
        hdu = [Mock(header=header)]
        self.assertEqual(util.frame_label(hdu), 2779494)

        # calexp-style: no EXPNUM, use FILENAME with trailing o removed
        header_calexp = Header([('FILENAME', '2779465o'), ('EXPID', 100060740)])
        hdu2 = [Mock(header=header_calexp)]
        self.assertEqual(util.frame_label(hdu2), 2779465)

        header_empty = Header([])
        hdu3 = [Mock(header=header_empty)]
        self.assertEqual(util.frame_label(hdu3, fallback_name='DIFFEXP-1111111-01.fits'), 1111111)

    def test_night_from_path(self):
        self.assertEqual(
            util.night_from_path('/data/AS1 July/2022-08-26/ccd00/DIFFEXP-2779494-00.fits'),
            '2022-08-26',
        )
        self.assertIsNone(util.night_from_path('/data/2779494/DIFFEXP-00.fits'))

    def test_stack_night_from_warpd(self):
        header = Header([
            ('WARPD', '/arc/classy/AS1 July/2022-08-26/ccd00'),
            ('MJD-OBS', 59817.4833963),
            ('MJDEND', 59817.4868922),
        ])
        hdu = [Mock(header=header)]
        self.assertEqual(util.stack_night(hdu), '2022-08-26')

    def test_stack_night_from_exposure_time(self):
        header = Header([
            ('MJD-OBS', 59817.4833963),
            ('MJDEND', 59817.4868922),
        ])
        hdu = [Mock(header=header)]
        self.assertEqual(util.stack_night(hdu), '2022-08-26')

    def test_field_from_path(self):
        path = '/data/AS1 July/2022-08-26/ccd00/DIFFEXP-2779494-00.fits'
        self.assertEqual(util.field_from_path(path), 'AS1 July')
        self.assertEqual(util.sanitize_filename_component('AS1 July'), 'AS1_July')

    def test_source_name_from_orbit(self):
        obs = Mock(provisional_name='2020 AB123')
        orbit = Mock(observations=[obs])
        orbit.name = 'LASTOBS'
        self.assertEqual(util.source_name_from_orbit(orbit), '2020 AB123')

        orbit_no_obs_name = Mock(observations=[])
        orbit_no_obs_name.name = 'FALLBACK'
        self.assertEqual(util.source_name_from_orbit(orbit_no_obs_name), 'FALLBACK')

        orbit_empty = Mock(observations=[Mock(provisional_name='   ')])
        orbit_empty.name = 'FALLBACK'
        self.assertEqual(util.source_name_from_orbit(orbit_empty), 'FALLBACK')

        orbit_np = Mock(observations=numpy.array([obs]))
        orbit_np.name = 'LASTOBS'
        self.assertEqual(util.source_name_from_orbit(orbit_np), '2020 AB123')

    def test_stack_output_filename(self):
        name = util.stack_output_filename(
            'AS1 July', '2022-08-26', 0, 1, 0, 'WEIGHTED_MEDIAN', 2.0, 45.0,
            source_name='2020 AB123',
        )
        self.assertEqual(
            name,
            'STACK-2020_AB123-AS1_July-2022-08-26-00-01-00-WEIGHTED_MEDIAN+02.00-+45.00.fits',
        )
        name_no_source = util.stack_output_filename(
            '03148', '2022-08-26', 12, 3, 1, 'MEDIAN', 1.5, -8.0,
        )
        self.assertEqual(
            name_no_source,
            'STACK-03148-2022-08-26-12-03-01-MEDIAN+01.50--08.00.fits',
        )

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
