#!/usr/bin/env python3

"""chromcorrect.py

To fit/correct 2d chromatograms

TODO: Rename to chrom-compose?
"""

import math
from subprocess import run
import numpy as np
import csv

import sys
import argparse

# from scipy import minimize


def csvWriter(filename, x, y):
    with open(filename, 'w') as f:
        writer = csv.writer(f)
        writer.writerows(zip(x, y))

def readfile(data_path, columns=[0,1], header=False, xticksColumn=0):
    """ Read x-y CSV-style files
    """
    if ':' in data_path:
        run(['scp', '-rC', data_path, '/tmp/plotting.csv'])
        data_path = '/tmp/plotting.csv'

    x = []
    y = []
    xticks = []
    # columns = [0, 1]
    delimiter = ' '
    with open(data_path, newline='') as csvfile:
        if ',' in csvfile.readline():
            delimiter = ','
    with open(data_path, newline='') as infile:
        # data = list(csv.reader(infile))
        if header:
            print(infile.readline())
        for line in infile:
            data_line = line.strip().split(delimiter)
            data_line = list(filter(None, data_line))
            if (data_line != []):
                if len(data_line) == 1:
                    y.append(float(data_line[0]))
                else:
                    x.append(float(data_line[columns[0]]))
                    y.append(float(data_line[columns[1]]))
                # if columns[0] != -1:
                #     x.append(float(data_line[columns[0]]))
                # y.append(float(data_line[columns[1]]))
                if xticksColumn is not None: 
                    xticks.append(float(data_line[xticksColumn]))

    return x, y, xticks


def getCentersAndAreas(col_radius:float, nrad:int, radial_disc_type='EQUIDISTANT'): 
    nRegions = nrad
    nShells = nRegions + 1 #Including r = 0
    rShells = []

    ## NOTE: move to current unit?
    # R = params.col_radius
    R = col_radius

    if radial_disc_type == 'EQUIVOLUME':
        for n in range(nShells):
            rShells.append(R * math.sqrt(n/nRegions))
    elif radial_disc_type == 'EQUIDISTANT':
        for n in range(nShells):
            rShells.append(R * (n/nRegions))

    centers = [ (rShells[i+1] + rShells[i])/2 for i in range(nRegions) ]
    areas = [ math.pi *  (rShells[i+1]**2 - rShells[i]**2) for i in range(nRegions) ]

    return centers, areas

def split_to_full(splitchromas, flowrates):

    denominator = sum(flowrates)
    numerator = np.sum(
            list(map(lambda x: [ y*x[1] for y in x[0]]  , zip(splitchromas, flowrates))),
            axis = 0)

    final = [ x/denominator for x in numerator]

    # times = readfile(files[0])[0]
    # csvWriter('composed_chromatogram.csv', times, final)
    return final

def evaluator(splitchromas, flowrates, concs): 
    denominator = sum(flowrates)
    numerator = np.sum(
            list(map(lambda x: [ y*x[1] for y in x[0]]  , zip(splitchromas, flowrates))),
            axis = 0)

    final = [ x/denominator for x in numerator ]

    return []

def main():

    ap = argparse.ArgumentParser()
    ap.add_argument('chromatograms', nargs='*', help='<time> <concentration> chromatogram csv files')
    ap.add_argument('-f', '--flowrates', help='<centers> <flowrates> csv file')

    args = ap.parse_args()

    files = args.chromatograms
    nrad = len(files)

    flowrates = readfile(args.flowrates)[1]
    assert len(flowrates) == nrad

    splitchromas = list(map(lambda f: readfile(f)[1], files))

    t,c,_ = readfile(files[0])

    concs = [ x * sum(flowrates) for x in c] 
    
    # splitchromas = [ [0] * len(c)  ] * nrad

    composed = split_to_full(splitchromas, flowrates)

    csvWriter('composed_chromatogram.csv', t, composed)


if __name__ == "__main__": 
    main()
