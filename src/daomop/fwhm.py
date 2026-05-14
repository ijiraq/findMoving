from photutils.psf import fit_fwhm
from astropy.nddata.utils import Cutout2D
import numpy as np

def measure_fwhm(data, position, size=21):
    cutout = Cutout2D(data, position, size=size)
    y_center, x_center = cutout.data.shape[0] / 2, cutout.data.shape[1] / 2
    fwhm = fit_fwhm(cutout.data, xypos=(x_center, y_center))
    if len(fwhm) != 1:
        raise ValueError("FWHM fitting did not return a single value.")
    return fwhm[0]


def aperture_correction(fwhm, aperture_radius):
    """
    Compute the aperture correction using a 4-parameter exponential model.

    Parameters:
    - fwhm (float): Full Width at Half Maximum of the PSF.
    - aperture_radius (float): Radius of the aperture used.

    Returns:
    - float: Aperture correction value.
    """
    x = 2 * aperture_radius / fwhm

    # Best-fit parameters from the model
    a = 2.22152
    b = 1.08226
    c = 3.81088
    d = 2.82771

    # 4-parameter exponential model
    return a * np.exp(-b * x) + c * np.exp(-d * x)