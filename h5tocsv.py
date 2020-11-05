#!/bin/env python3

"""
Extracts csv files from the given h5 files.

@usage: `./h5tocsv.py <h5 files list> -e <solution path lists>`
        where solution paths provided are appended to sim.root.output.solution

if the given paths do not exist, or are empty, no csv is generated.
"""

from cadet import Cadet
import csv
import argparse

cadetpath = "/home/jayghoshter/local/bin/cadet-cli"
Cadet.cadet_path = cadetpath


ap = argparse.ArgumentParser()
ap.add_argument("files", nargs='*', help="h5 files")
ap.add_argument("-e", "--extract", nargs='*')
args = vars(ap.parse_args())


for f in args['files']:
    sim = Cadet()
    sim.filename = f
    sim.load()
    for varstring in args['extract']:
        var = varstring.split('/')
        x = sim.root.output.solution.solution_times
        y = sim.root.output.solution
        for item in var:
            y = y[item]
        if y == [] or y == {}:
            continue
        with open(f + '.' + varstring.replace('/', '_') +'.csv', 'w') as f:
            writer = csv.writer(f, delimiter=',')
            writer.writerows(zip(x, y))
