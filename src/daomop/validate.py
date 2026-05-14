import pyds9
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.io import fits

class ImageDisplay:
    
    def __init__(self, image_path):
        self.image_path = image_path
        self.ds9 = pyds9.DS9()
        self.ds9.set('file {}'.format(image_path))
        self.wcs = WCS(fits.open(image_path)[0].header)

    def mark_regions(self, coordinates, radius=5):
        """
        Mark regions on the DS9 display based on input list of RA/DEC coordinates using circles.
        
        Parameters:
        coordinates (list of tuples): List of (RA, DEC) coordinates to mark.
        radius (int): Radius of the circles to mark.
        """
        for ra, dec in coordinates:
            sky_coord = SkyCoord(ra, dec, unit='deg', frame='icrs')
            x, y = sky_coord.to_pixel(self.wcs)
            self.ds9.set('regions', 'image; circle({}, {}, {})'.format(x, y, radius))

    def get_imexamine_locations(self):
        """
        Returns the x/y locations returned by DS9 imexamine.
        
        Returns:
        list of tuples: List of (x, y) coordinates.
        """
        self.ds9.set('imexam')
        coords = []
        while True:
            result = self.ds9.get('iexam key coordinate')
            if result == 'q':
                break
            x, y = map(float, result.split())
            coords.append((x, y))
        return coords

# Example usage:
# display = ImageDisplay('path/to/image.fits')
# display.mark_regions([(10.684, 41.269), (83.822, -5.391)], radius=10)
# locations = display.get_imexamine_locations()
# print(locations)
