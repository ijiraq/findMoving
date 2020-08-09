import argparse

base_parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                      fromfile_prefix_chars='@',
                                      add_help=False)
base_parser.add_argument('basedir', help="Root directory of LSST pipelined data")
base_parser.add_argument('--pointing', help="Which pointing to process, eg. 03071", default=None)
base_parser.add_argument('--field', help='Which FIELD to process', choices=['NHF1', 'NHF2'], default=None)
base_parser.add_argument('--patch', help="sky patch to process (e.g. 0,0)", nargs=1)
base_parser.add_argument('--rerun', help="rerun directory containing the warped difference images.", nargs=1)
base_parser.add_argument('--filter', help="Filter to stack", default="HSC-R2")
base_parser.add_argument('--visit', help='visit to process.', default=None)
base_parser.add_argument('--ccd', help="Which CCD to stack?", type=int, default=None)
base_parser.add_argument('--log-level', help="What level to log at? (ERROR, INFO, DEBUG)", default="ERROR",
                         choices=['INFO', 'ERROR', 'DEBUG'])


def parse_rerun(rerun):
    reruns = rerun[0].split(":")
    if len(reruns) > 2:
        raise ValueError("Don't know what to do with more then 2 rerun directories.")

    if len(reruns) == 1:
        input_rerun = output_rerun = reruns[0]
    else:
        input_rerun = reruns[0]
        output_rerun = reruns[1]
    return input_rerun, output_rerun
