import argparse
import glob
import logging
import os

import numpy
from astropy import units
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata.utils import Cutout2D, PartialOverlapError
from astropy.table import Table, vstack
from astropy.wcs import wcs

PIX_CUTOUT_SIZE = 64


def load_plant_list(filename):
    """
    Load a plantList file as formatted from the planting pipeline for the NH HSC search.

    :param filename:  name of the .plantList file
    :return: table of values
    :rtype Table
    """

    column_names = ['index', 'ra', 'dec', 'x', 'y', 'rate', 'angle', 'rate_ra', 'rate_dec',
                    'mag', 'psf_amp']
    column_units = {'index': None,
                    'ra': units.degree,
                    'dec': units.degree,
                    'x': None,
                    'y': None,
                    'rate': units.arcsecond / units.hour,
                    'rate_ra': units.arcsecond / units.hour,
                    'rate_dec': units.arcsecond / units.hour,
                    'mag': None,
                    'psf_amp': None}
    this_table = Table.read(filename, names=column_names, format='ascii')
    for name in column_units:
        if column_units[name] is not None:
            this_table[name] = this_table[name] * column_units[name]
    return this_table


def cut(pairs, plant_list, random=True, size=PIX_CUTOUT_SIZE, num_samples=1000, extno=1):
    """
    retrieve image sections from two images.

    :param pairs: list of pairs of FITS filenames containing the images of interest.
    :param plant_list: list of sources that were added to these images.
    :param random: should the centre of the cutouts be random?
    :param size: dimensions of the cutout of 2D image (square cutout)
    :param num_samples: number of samples to return, if random (otherwise does full grid).
    :param extno: which extension in the FITS file contains the image data.
    :return: Three dictionaries, first contains cutouts with sources, then source locations, then cutouts without source
    :rtype list, list, list
    """

    # for each pair of images in the list of pairs make cutouts, either on a grid or at random
    # and record a data array that packs the two cutout sections beside each other.
    # Only do this if one or zero artificial sources are in that part of the image, if two
    # artificial sources appear, skip this cutout.
    # Also record the X, Y, XO, YO, MAG of the planted source for each cutout.
    # X/Y are in the original reference frame and XO/YO in the cutout reference frame.

    #
    source_cutouts = []
    blank_cutouts = []
    source_cutout_targets = []

    for pair in pairs:

        logging.info("Creating {} cutouts for image pair {}-{}".format(num_samples, pair[0], pair[1]))
        images = []
        pair_source_cutouts = []
        pair_blank_cutouts = []
        pair_targets = []

        # open the two FITS images that are the pair.
        for filename in pair:
            images.append(fits.open(filename)[extno])

        # Require that the images are the same dimensions.
        if images[0].data.shape != images[1].data.shape:
            raise IOError("Images must be the same dimension and registered.")

        # Make a list of X/Y coordinates to for the centres of cutouts to make
        # from the images in this pair.
        if random:
            xx = numpy.random.randint(0, images[0].header['NAXIS1'], num_samples)
            yy = numpy.random.randint(0, images[0].header['NAXIS2'], num_samples)
        else:
            xx = numpy.arange(size//2, images[0].header['NAXIS1'], size)
            yy = numpy.arange(size//2, images[0].header['NAXIS2'], size)

        # loop over the mess of coordinates.
        for x in xx:
            for y in yy:
                try:
                    image_cutouts = []
                    image_targets = []
                    use_this_cutout = False
                    for image in images:
                        image_wcs = wcs.WCS(image.header)
                        image_cutout = Cutout2D(image.data, (x, y), size, wcs=image_wcs, mode='strict')
                        # determine which planted sources are in this cutout.
                        image_fakes = plant_list[image_cutout.wcs.footprint_contains(plant_list['skycoord'])]
                        if len(image_fakes) > 1:
                            # only use image cutouts with 0 or 1 source.
                            break
                        image_cutouts.append(image_cutout.data)
                        target_xo = target_yo = target_x = target_y = target_mag = -1
                        use_this_cutout = True
                        if len(image_fakes) == 1:
                            cutout_coo = image_cutout.wcs.all_world2pix(((image_fakes['ra'][0],
                                                                        image_fakes['dec'][0]),),
                                                                        1)
                            target_xo = cutout_coo[0][0]
                            target_yo = cutout_coo[0][1]

                            initial_coo = image_wcs.all_world2pix(((image_fakes['ra'][0],
                                                                   image_fakes['dec'][0]),),
                                                                  1)
                            target_x = initial_coo[0][0]
                            target_y = initial_coo[0][1]
                            target_mag = image_fakes['mag'][0]
                        image_targets.append([target_x, target_y, target_xo, target_yo, target_mag])
                    if use_this_cutout:
                        if numpy.sum(image_targets[0]) > 0 or numpy.sum(image_targets[1]) > 0:
                            pair_source_cutouts.extend([image_cutouts])
                            pair_targets.extend([image_targets])
                        else:
                            pair_blank_cutouts.extend([image_cutouts])
                except PartialOverlapError as ex:
                    logging.debug("Cutout at {},{} not fully inside image, skipping".format(x, y))
                    # skip if only a partially overlapping cutout.
                    break

        # append data to .npy file for this pair.
        # when reading these npy files we loop over them until we get a
        # read exception.

        source_cutouts.append(pair_source_cutouts)
        source_cutout_targets.append(pair_targets)
        blank_cutouts.append(pair_blank_cutouts)
    return source_cutouts, source_cutout_targets, blank_cutouts


def build_image_pair_list(image_directory):
    """
    Build a list of images to make cutouts from and return list of pairs of images.
    :param image_directory:
    :return: list of image pairs
    :rtype list
    """
    image_filename_list = numpy.array(glob.glob(os.path.join(image_directory, '*.fits')))
    image_pairs = numpy.random.choice(image_filename_list,
                                      (len(image_filename_list)//2, 2),
                                      replace=False)
    return image_pairs


def build_table_of_planted_sources(plant_list_directory, pattern='*.plantList'):
    """
    Build a table containing all the sources in .plantList files in the given directory.

    This method calls load_plant_list and expects a plantList file in a specific format.

    :param plant_list_directory: directory containing plant list files.
    :param pattern: pattern to match against (defaults to *.plantList)
    :return: Table of sources from all the plantList files in plant_list_directory whose filenames match pattern
    :rtype Table
    """

    plant_filename_list = glob.glob(os.path.join(plant_list_directory, pattern))
    table_of_planted_sources = None
    for plant_filename in plant_filename_list:
        if table_of_planted_sources is None:
            table_of_planted_sources = load_plant_list(plant_filename)
        else:
            logging.debug("Reading artificial sources from {}".format(plant_filename))
            table_of_planted_sources = vstack((table_of_planted_sources, load_plant_list(plant_filename)))

    # create a skycoord column, for later use.
    table_of_planted_sources['skycoord'] = SkyCoord(table_of_planted_sources['ra'],
                                                    table_of_planted_sources['dec'])

    return table_of_planted_sources


def main():
    """"
    This is the main method which runs this module as a Command Line script.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--image-directory', default='./',
                        help='name of directory containing images to draw sections from.')
    parser.add_argument('--plant-list-directory', default='./',
                        help='name of directory where plantList files can be found.')
    parser.add_argument('--nsamples', type=int,
                        default=1000, help='How many cutouts to draw from the image.')
    parser.add_argument('--random', action='store_true', default=False,
                        help='Do cutouts as random locations, instead of grid')
    parser.add_argument('--dimension', type=int, default=PIX_CUTOUT_SIZE,
                        help='width/height of square cutouts to make.')
    parser.add_argument('--verbose', default=False, action='store_true')

    args = parser.parse_args()
    if args.verbose:
        logging.basicConfig(level=logging.INFO)
    else:
        logging.basicConfig(level=logging.ERROR)

    # Get the list of pairs of images in image_directory
    image_pairs = build_image_pair_list(args.image_directory)
    logging.info("Created a list of {} image pairs.".format(len(image_pairs)))
    table_of_planted_sources = build_table_of_planted_sources(args.plant_list_directory)
    logging.info("Read in {} artificial sources.".format(len(table_of_planted_sources)))
    logging.info("Extracting {} cutouts of size {}X{} from those image pairs".format(args.nsamples,
                                                                                     args.dimension,
                                                                                     args.dimension))
    if args.random:
        logging.info('Grid will be random')
    else:
        logging.info('Grid will be regularized')

    source_cutouts, source_targets, blank_cutouts = cut(image_pairs, table_of_planted_sources,
                                                        size=args.dimension,
                                                        random=args.random,
                                                        num_samples=args.nsamples)
    # print the location of the planted source in the first image of the first pair of images
    # along with the shape of the data array of the cutout from the first image that contains that source
    # Indexing is 'Image_Pair', 'Cutout Set', 'image1, image2'
    print(source_targets[0][0][0][2:4], source_cutouts[0][0][0].shape)
    # print the location of the planted source in the second image of the first pair of images
    # along with the shape of the data array of the cutout from the second image that contains that source
    print(source_targets[0][0][1][2:4], source_cutouts[0][0][1].shape)
    # print out the shape of the first blank piece of sky (no artificial TNO) from the first image in the first pair.
    print(blank_cutouts[0][0][0].shape)
    # print out the shape of the first blank piece of sky (no artificial TNO) from the first image in the first pair.
    print(blank_cutouts[0][0][1].shape)


if __name__ == '__main__':
    main()