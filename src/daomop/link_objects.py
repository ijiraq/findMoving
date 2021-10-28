"""
Given a target .mpc file an a setof possible additoinal measures of target find the best matching link..

Links are consider valid if:
    - candidate link file contains mpc lines.
    - adding those observations results in a bound orbit with all ra/dec residuals less than 0.2 arc-seconds
"""
import mp_ephem
from astropy import units
import logging
import argparse
import os

ObsRecord = mp_ephem.ObsRecord
BKOrbit = mp_ephem.BKOrbit


def main():
    """
    This is the function that runs if you use this module as a console script

    :return: Number of candidates that provided a link.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('target', help="the file containing MPC lines of the thing you want to find another observations of")
    parser.add_argument('candidates', nargs='*', help="files that contain MPC lines of detected sources that might link to target")
    parser.add_argument('--log-level', default='INFO', choices=['DEBUG', 'INFO', 'ERROR'])
    args = parser.parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level))
    return link(args.target, args.candidates)


def link(target, candidates):
    """

    :param target: file containing MPC lines of thing to try and link to.
    :param candidates: list of files containing MPC lines of sources that might be matches.
    :return: count of how many candidates provided a link.
    """

    # Load the observations of the target into a list.
    baseObservations = []
    nlinks = 0
    with open(target) as fobj:
        for line in fobj.readlines():
            baseObservations.append(ObsRecord.from_string(line))
    logging.debug(f"Retrieved {len(baseObservations)} for {target}")
    if not len(baseObservations) > 0:
        return -1

    for filename in candidates:
        logging.debug(f"Attempting to link {filename}")
        new_obs = []
        linkable_candidate = True
        with open(filename) as fobj:
            for line in fobj.readlines():
                obs = ObsRecord.from_string(line)
                if obs.null_observation:
                    logging.debug("Not keeping NULL observations.")
                    continue
                new_obs.append(obs)
        if not len(new_obs) > 0:
            continue

        trial_obs = baseObservations + new_obs
        try:
            orbit = BKOrbit(trial_obs)
            orbit.predict("2021-01-01T00:00:00")
            orbit.compute_residuals()
        except Exception as ex:
            logging.error(f"Linking attempting link with {filename} -> {ex}")
            continue

        ll = orbit

        # Check all the positions for goodness of fit
        for obs in orbit.observations:
            if obs.ra_residual > 0.2 or obs.dec_residual > 0.2:
                linkable_candidate = False

        if not linkable_candidate:
            logging.debug(f"Link with {filename} results in large residuals, rejecting")
            continue

        if orbit.a < 0 * units.au or orbit.inc > 90 * units.degree:
            logging.debug(f"Link with {filename} results in non-physical orbit, rejecting")
            continue

        if linkable_candidate:
            print(f"{target} + {filename}")
            print(ll.summarize())
            nlinks += 1
            with open(f"{os.path.splitext(target)[0]}_{os.path.splitext(filename)[0]}.mpc", 'w') as fobj:
                for obs in ll.observations:
                    fobj.write(obs.to_string()+'\n')

    return nlinks


if __name__ == '__main__':
    print(main())
