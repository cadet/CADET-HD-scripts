#!/bin/env python3

"""
Extracts csv files from the given h5 files.

@usage: `./h5tocsv.py <h5 files list> -e <solution path lists>`
        where solution paths provided are appended to sim.root.output.solution

if the given paths do not exist, or are empty, no csv is generated.
"""

import itertools

from cadet import Cadet
import csv
import argparse
import numpy as np

cadetpath = "/home/jayghoshter/local/bin/cadet-cli"
Cadet.cadet_path = cadetpath


ap = argparse.ArgumentParser()
ap.add_argument("-f", "--files", nargs='*', help="h5 files", required=True)
ap.add_argument("-e", "--extract", nargs='*')
ap.add_argument("-i", "--integrate", nargs='*')
args = vars(ap.parse_args())


for f in args['files']:
    sim = Cadet()
    sim.filename = f
    sim.load()

    if args['extract']:
        for varstring in args['extract']:
            var = varstring.split('/')
            x = sim.root.output.solution.solution_times
            y = sim.root.output.solution
            for item in var:
                y = y[item]
            if len(y) == 0:
                continue

            ## flatten list (in case of bulk solution)
            while isinstance(y[0], np.ndarray):
                flatten = itertools.chain.from_iterable
                y = list(flatten(y))

            with open(''.join([f, '.', varstring.replace('/', '_'), '.csv' ]), 'w') as outfile:
                writer = csv.writer(outfile, delimiter=',')
                writer.writerows(zip(x, y))

    if args['integrate']:
        """WARNING
            [NOTE] Case-specific code follows
        """
        for varstring in args['integrate']:
            var = varstring.split('/')
            x = sim.root.output.solution.solution_times
            y = sim.root.output.solution
            for item in var:
                y = y[item]
            if len(y) == 0:
                continue

            # ndims = len(y.shape)
            # for i in range(ndims-2,0,-1):
            unit = varstring.split('/')[0]
            print(unit)

            if sim.root.input.model[unit].discretization.par_disc_type == b'EQUIDISTANT_PAR':
                par_radius = sim.root.input.model[unit].par_radius
                npar = sim.root.input.model[unit].discretization.npar
                ## WARNING: Probably wrong?
                dx_particle = par_radius/(npar-1)
                print("par radius:", par_radius)
                print("npar: ", npar)
                print("dx particle:", dx_particle)
            else:
                print(sim.root.input.model[unit].discretization.par_disc_type)
                raise NotImplementedError

            dx_column = sim.root.input.model[unit].discretization.ncol

            print(varstring)
            print("Shape of Y: {shape}".format(shape=y.shape))
            print("Time | Axial | Particle | ...")

            y = np.trapz(y, dx=dx_particle, axis=2)
            # y = np.trapz(y, dx=dx_column, axis=1)

            par_radius = sim.root.input.model[unit].par_radius

            ## TODO: Use numpy.squeeze to flatten array axis with length = 1

            y = [ g[0][0] / (par_radius) for g in y ]
            # y = [ g[0][0] for g in y ]

            print(y)

            with open(f + '.integrated.' + varstring.replace('/', '_') +'.csv', 'w') as fd:
                writer = csv.writer(fd, delimiter=',')
                writer.writerows(zip(x, y))
