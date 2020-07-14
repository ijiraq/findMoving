import argparse
import glob
import logging
import os
import re
import sqlite3
from itertools import combinations

import numpy
from astropy import units
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata.utils import Cutout2D, PartialOverlapError
from astropy.table import Table
from astropy.visualization import ImageNormalize, ZScaleInterval
from astropy.wcs import wcs
from matplotlib import pyplot

PIX_CUTOUT_SIZE = 64
VISIT_IN_PLANT_LIST_FILENAME_RE = re.compile(r'0([0-9]{6})')
PLANT_LIST_COLUMNS = ['index', 'ra', 'dec', 'x', 'y', 'rate', 'angle', 'rate_ra', 'rate_dec', 'mag', 'psf_amp']
PLANT_LIST_UNITS = {'index': None,
                    'ra': units.degree,
                    'dec': units.degree,
                    'x': units.pix,
                    'y': units.pix,
                    'rate': units.arcsecond / units.hour,
                    'angle': units.degree,
                    'rate_ra': units.arcsecond / units.hour,
                    'rate_dec': units.arcsecond / units.hour,
                    'mag': units.mag,
                    'psf_amp': units.adu}
PLANT_LIST_db = {'visit': 'INTEGER',
                 'index': 'INTEGER',
                 'ra': 'REAL',
                 'dec': 'REAL',
                 'x': 'REAL',
                 'y': 'REAL',
                 'rate': 'REAL',
                 'angle': 'REAL',
                 'rate_ra': 'REAL',
                 'rate_dec': 'REAL',
                 'mag': 'REAL',
                 'psf_amp': 'REAL'
                 }


def insert_plant_list_into_database(this_table, plant_list_db='plant_list.db'):
    """

    :param this_table: the plantList format astropy.table.Table to load into the SQL DB.
    :param plant_list_db: the sqlite3 database to load them into.
    :return: None
    """
    # Compare the first line of the file I've been given to be sure it matches what is expected.
    if not os.access(plant_list_db, os.R_OK):
            init_db(plant_list_db)
    cast={'REAL': float, 'INTEGER': int }
    with sqlite3.connect(plant_list_db) as db:
        visit = this_table['visit'][0]
        cursor = db.cursor()
        # Check if we've already added this visit
        cursor.execute(f'SELECT count(*) FROM fakes WHERE visit={visit}')
        row = cursor.fetchone()
        if row[0] > 0:
            logging.warning(f'Already ingested {visit} in to {plant_list_db}')

        sql = 'REPLACE INTO fakes(' + ",".join([f'`{name}`' for name in this_table.colnames]) + ')'
        sql += ' VALUES(' + ','.join(['?']*len(this_table.colnames)) + ')'
        logging.debug(f'INSERTING using:\n{sql}')
        for row in this_table:
            values = [cast[PLANT_LIST_db[name]](row[name]) for name in this_table.colnames]
            logging.debug(f'values:\n{values}')
            db.execute(sql, values)
        db.commit()
    return


def get_visit_plant_list(visit, plant_list_db):
    """
    Retrieve all the sources associated with the given visit into a astropy.table.Table

    :param plant_list_db: sqlite3 database (likely created with init_db method in this module)
    :return: Table of planted sources in plantList format.
    """
    with sqlite3.connect(plant_list_db) as db:
        cursor = db.cursor()
        sql = f'SELECT * FROM `fakes` WHERE `visit`={visit}'
        result = numpy.array(cursor.execute(sql).fetchall())
        colnames = [x[0] for x in cursor.description]
    # logging.debug(f'Loaded columns: {colnames}')
    # logging.debug(f'Table shape: {result.shape}')
    this_table = Table(data=result, names=colnames)
    this_table['skycoord'] = SkyCoord(this_table['ra'], this_table['dec'], unit=('degree', 'degree'))
    return this_table


def load_plantlist(filename):
    """
    Load a plantList file as formatted from the planting pipeline for the NH HSC search.

    These files were created by Wes and have a funny format.

    :param filename:  name of the .plantList file
    :return: table of values
    :rtype Table
    """
    logging.debug(f'Loading plantList from {filename}')
    plantlist_first_line = '#index ra dec x y rate ("/hr) angle (deg) rate_ra rate_dec mag psf_amp'
    first_line_in_file = open(filename, 'r').readline().strip()
    assert first_line_in_file == plantlist_first_line, \
        'expected:{}\ngot:{}\n is not a plantList formatted file.'.format(plantlist_first_line,
                                                                          first_line_in_file,
                                                                          filename)
    column_names = PLANT_LIST_COLUMNS
    column_units = PLANT_LIST_UNITS
    match = VISIT_IN_PLANT_LIST_FILENAME_RE.search(filename)

    if match is None:
        raise ValueError(f'Filename {filename} does not match expected pattern.')
    visit = int(match.group(1))

    this_table = Table.read(filename, names=column_names, format='ascii')
    for name in column_units:
        if column_units[name] is not None:
            this_table[name] = this_table[name] * column_units[name]

    this_table['visit'] = visit

    return this_table


def cut(pairs, full_plant_list, random=True, size=PIX_CUTOUT_SIZE, num_samples=1000, extno=1):
    """
    retrieve image sections from image sets whose files names are given in pairs in the pairs list.

    :param pairs: list of pairs of FITS filenames containing the images of interest.
    :param full_plant_list: list of sources that were added to these images.
    :param random: should the centre of the cutouts be random?
    :param size: dimensions of the cutout of 2D image (square cutout)
    :param num_samples: number of samples to return, if random (otherwise does full grid).
    :param extno: which extension in the FITS file contains the image data.
    :return: Three lists, first contains cutouts with sources, then source locations, then cutouts without a source
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

        logging.info("Creating {} cutouts for image pair {}".format(num_samples, pair))
        images = {}

        # open the two FITS images that are the pair.
        shape = None
        for filename in pair:
            visit = re.search(r'([0-9]{6})', filename).group(1)
            images[visit] = fits.open(filename)[extno]
            if shape is None:
                shape = images[visit].shape
            assert shape == images[visit].shape, "All images must be the same dimension and registered."
            plant_list = get_visit_plant_list(visit, full_plant_list)
            logging.info("{} sources planted into visit {}".format(len(plant_list), visit))

        # Make a list of X/Y coordinates to for the centres of cutouts to make
        # from the images in this pair.
        visit = list(images.keys())[0]
        if random:
            xx = numpy.random.randint(size // 2, images[visit].header['NAXIS1'] - size // 2, num_samples * 2)
            yy = numpy.random.randint(size // 2, images[visit].header['NAXIS2'] - size // 2, num_samples * 2)
        else:
            xx = numpy.arange(size // 2, images[visit].header['NAXIS1'], size)
            yy = numpy.arange(size // 2, images[visit].header['NAXIS2'], size)
            num_samples = len(xx)
        # Make an array that contains the coordinate pairs centred on x[o],y[o].
        xy = numpy.vstack((xx, yy)).T
        num_cutouts = 0
        num_cutouts_with_source = 0
        # For the 'random' set we pad out the number of samples to allow for nan values being skipped.
        # loop over the mess of coordinates.
        logging.info(f'Looping over {len(xy)} coordinate pairs to make cutouts.')
        for p in xy:
            try:
                image_cutouts = []
                image_targets = []
                use_this_cutout = False
                for visit in images:
                    image = images[visit]
                    plant_list = get_visit_plant_list(visit, full_plant_list)
                    image_wcs = wcs.WCS(image.header)
                    image_cutout = Cutout2D(image.data, p, size, wcs=image_wcs, mode='strict')
                    image_cutout.data = numpy.nan_to_num(image_cutout.data, nan=-100)
                    if numpy.isnan(image_cutout.data).any():
                        raise ValueError("Cutout of {} at {} has NaN".format(p, visit))
                    # determine which planted sources are in this cutout.
                    image_fakes = plant_list[image_cutout.wcs.footprint_contains(plant_list['skycoord'])]
                    if len(image_fakes) > 1:
                        raise ValueError("Cutout at {} of {} has too many ({}) fakes.".format(p, visit,
                                                                                              len(image_fakes)))
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
                    image_targets.append(numpy.array([target_x, target_y, target_xo, target_yo, target_mag]))
                if use_this_cutout and len(image_targets) == len(pair) and len(image_cutouts) == len(pair):
                    if numpy.sum(image_targets[0]) > 0 or numpy.sum(image_targets[1]) > 0:
                        source_cutouts.append(numpy.array(image_cutouts))
                        source_cutout_targets.append(numpy.array(image_targets))
                        num_cutouts_with_source += 1
                    else:
                        blank_cutouts.append(numpy.array(image_cutouts))
                        logging.debug("blank shape: {}".format(blank_cutouts[-1].shape))
                    num_cutouts += 1
            except PartialOverlapError as ex:
                # skip this cutout if it is only a partially overlapping cutout.
                logging.debug(str(ex))
                logging.debug("Cutout at {} not fully inside image, skipping".format(p))
            except ValueError as io:
                logging.debug(str(io))
            if num_cutouts > num_samples:
                break
        logging.info(f'Extracted {num_cutouts} from {pair} with {num_cutouts_with_source} '
                     f'containing an artificial source.')

    # restructure the data so that we loose the 'pair' breakdown.
    source_cutouts = numpy.array(source_cutouts)
    source_cutout_targets = numpy.array(source_cutout_targets)
    source_cutout_targets = numpy.array(source_cutout_targets)
    blank_cutouts = numpy.array(blank_cutouts)
    logging.debug("Sending back data with the following shapes: {}  {}  {}".format(source_cutouts.shape,
                                                                                   source_cutout_targets.shape,
                                                                                   blank_cutouts.shape))
    return source_cutouts, source_cutout_targets, blank_cutouts


def build_image_pair_list(image_directory, num_pairs=None, random=False, fraction=1.0, num_per_pair=2):
    """
    Build a list of images to make cutouts from and return list of pairs of images.
    :param image_directory:
    :param num_pairs: How many image pairs to build? None ==> All possible.
    :param num_per_pair: How many images in a 'pair' (normally 2)
    :param random: should the pairs be random or an iteration over all combinations.
    :param fraction: if random, what fraction of a complete sample should be provided (there will be repeats)?
    :return: list of image pairs
    """
    image_filename_list = numpy.array(glob.glob(os.path.join(image_directory, '*.fits')))
    logging.debug(f'Got list {image_filename_list} of file in directory {image_directory}')
    if random:
        # Make the pairs via random sampling
        image_pairs = numpy.random.choice(image_filename_list,
                                          int(fraction * len(image_filename_list) * (len(image_filename_list) - 1),
                                              num_per_pair),
                                          replace=True)
    else:
        image_pairs = combinations(image_filename_list, num_per_pair)
    if num_pairs is not None:
        image_pairs = [x for x in image_pairs][0:num_pairs]
    logging.debug(f'Sending back {len(image_pairs)} pairs of images.')
    return image_pairs


def build_table_of_planted_sources(plant_list_directory, pattern='*.plantList', plant_list_db=None, reload=False):
    """
    Build a table containing all the sources in .plantList files in the given directory.

    This method calls load_plantlist and expects a plantList file in a specific format.

    :param plant_list_directory: directory containing plant list files.
    :param pattern: pattern to match against (defaults to *.plantList)
    :return: Table of sources from all the plantList files in plant_list_directory whose filenames match pattern
    :rtype Table
    """
    if plant_list_db is None:
        plant_list_db = os.path.join(plant_list_directory, 'plant_list.db')
    if os.access(plant_list_db, os.R_OK) and not reload:
        return plant_list_db

    plant_filename_list = glob.glob(os.path.join(plant_list_directory, pattern))
    for plant_filename in plant_filename_list:
        try:
            this_table = load_plantlist(plant_filename)
            insert_plant_list_into_database(this_table,
                                            plant_list_db=plant_list_db)
        except Exception as ex:
            logging.error(str(ex))
            logging.warning("Skipping {}".format(plant_filename))

    # create a skycoord column, for later use.
    return plant_list_db


def plot(source_cutout, source_cutout_target):
    """
    Make a sky plot of the two image chunks in the source cutout and display for user.

    :param source_cutout: a single source cutout entry as returned from data_model.cut
    :param source_cutout_target: a single source cutout target entry as returned from data_model.cut
    :return: None
    """

    fig = pyplot.figure(figsize=(11, 11))

    for i in range(len(source_cutout)):
        fig.add_subplot(1, 2, i + 1)
        try:
            norm = ImageNormalize(source_cutout[i], interval=ZScaleInterval())
        except Exception as ex:
            logging.debug(str(ex))
            norm = None

        pyplot.imshow(source_cutout[i], origin='lower', norm=norm)
        # Offset the x/y location by 1 pixel due to locations being FORTRAN/FITS and plot being python/numpy.
        pyplot.scatter(source_cutout_target[i][2] - 1,
                       source_cutout_target[i][3] - 1,
                       s=500, facecolors='none', edgecolors='r')
        pyplot.title("planted source mag: {:4.1f}".format(source_cutout_target[i][4]))

    pyplot.show()


def main():
    """"
    This is the main method which runs this module as a Command Line script.
    """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--image-directory', default='./',
                        help='name of directory containing images to draw sections from.')
    parser.add_argument('--reload-plant-list', action='store_true',
                        help="Rebuild the plant_list database if it exists?")
    parser.add_argument('--plant-list-directory', default='./',
                        help='name of directory where plantList files can be found.')
    parser.add_argument('--nsamples', type=int,
                        default=1000, help='How many cutouts to draw from the image.')
    parser.add_argument('--npairs', type=int,
                        default=19, help='How many pairs of images do you want to build nsample cutouts of.')
    parser.add_argument('--num-per-pair', default=2, help='how many images in a pair, normally 2')
    parser.add_argument('--random', action='store_true', default=False,
                        help='Do cutouts as random locations, instead of grid')
    parser.add_argument('--dimension', type=int, default=PIX_CUTOUT_SIZE,
                        help='width/height of square cutouts to make.')
    parser.add_argument('--verbose', default=False, action='store_true')
    parser.add_argument('--num-to-plot', default=15, type=int,
                        help="Should we plot some example cutouts? How many? Set to 0 if you don't want any.")
    parser.add_argument('--debug', default=False, action='store_true')

    args = parser.parse_args()

    # Set the logging level.
    level = logging.ERROR
    if args.verbose:
        level = logging.INFO
    if args.debug:
        level = logging.DEBUG
    logging.basicConfig(level=level)

    # Get the list of pairs of images in image_dir
    image_pairs = build_image_pair_list(args.image_directory, num_pairs=args.npairs, num_per_pair=args.num_per_pair)
    logging.info("Created a list of image pairs.".format(image_pairs))
    plant_list_db = build_table_of_planted_sources(args.plant_list_directory, reload=args.reload_plant_list)
    # logging.info("Read in {} artificial sources.".format(len(table_of_planted_sources)))
    logging.info("Extracting {} cutouts of size {}X{} from those image pairs".format(args.nsamples,
                                                                                     args.dimension,
                                                                                     args.dimension))
    if args.random:
        logging.info('Grid will be random')
    else:
        logging.info('Grid will be regularized')

    source_cutouts, source_targets, blank_cutouts = cut(image_pairs, plant_list_db,
                                                        size=args.dimension,
                                                        random=args.random,
                                                        num_samples=args.nsamples)
    #
    # print the location of the planted source in the first image of the first pair of images
    # along with the shape of the data array of the cutout from the first image that contains that source
    # Indexing is 'Image_Pair', 'Cutout Set', 'image1, image2'
    #
    logging.info("source target shape {}.".format(source_targets.shape))
    logging.info("source cutouts shape {}".format(source_cutouts.shape))
    logging.info("blank cutouts shape {}.".format(blank_cutouts.shape))

    #
    # for the first image pair (index 0 on source_targets has all the planted source parameters for the
    # first pair of images). sort the cutout set (second index, given as : below selects all the cutouts for this pair)
    # on the magnitude of the planted source in first image in the pair (the next '0' in indexing says we are dealing
    # with the first image in the pairs), lower value of magnitude (stored in column 4) are brighter.
    #
    if args.num_to_plot > 0:
        brightest_targets = numpy.argsort(source_targets[:, 0, 4])
        # Plot cutout images of the 15 brightest sources.
        for i in range(min(len(brightest_targets), args.num_to_plot)):
            plot(source_cutouts[brightest_targets[i]], source_targets[brightest_targets[i]])


def init_db(dbfilename):
    """
    Create an SQLite database to add plantList files into.

    :param dbfilename:
    :return:
    """
    try:
        conn = sqlite3.connect(dbfilename)
        c = conn.cursor()
        # Create table, append to the list of columns in plantList a visit column.
        column_defs = ",".join([f'`{col}` {PLANT_LIST_db[col]}' for col in PLANT_LIST_db])
        sql = f'CREATE TABLE IF NOT EXISTS `fakes`({column_defs}, UNIQUE(`visit`,`index`))'
        logging.debug(f'Creating sql TABLE with: \n{sql}')
        c.execute(sql)
        conn.commit()
        conn.close()
    except Exception as ex:
        logging.error(str(ex))
        raise ex


if __name__ == '__main__':
    main()
