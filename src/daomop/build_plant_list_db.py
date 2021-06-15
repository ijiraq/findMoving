from daomop import data_model
import logging
import argparse


def main():
    parser = argparse.ArgumentParser(description='scan through all the plantList files in a directory and'
                                                 'create an sqlite3 DB containing those fake sources.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('plant_list_dir', help='Directory containing the .plantList files.')
    parser.add_argument('--plant-list-db', help='Name of the database to create.', default='plant_list.db')
    parser.add_argument('--log-level', default='INFO', choices=['DEBUG', 'INFO', 'ERROR'])
    parser.add_argument('--reload', action='store_true')
    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level))
    plant_list_directory = args.plant_list_dir
    data_model.build_table_of_planted_sources(plant_list_directory, plant_list_db=args.plant_list_db, 
                                              reload=args.reload)


if __name__ == '__main__':
    main()
