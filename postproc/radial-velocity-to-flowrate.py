#!/usr/bin/env python3

"""
Given an input of radial velocity profile, 
generate the corresponding flowrate profile.
"""

import csv
import argparse
from subprocess import run
import numpy as np
from scipy.interpolate import interp1d
import math

def extrapolate(x0, y0, x, kind='linear'): 
    f = interp1d(x0, y0, kind=kind, fill_value='extrapolate')
    return f(x)

def readfile(data_path, columns=[0,1], header=False):
    """ Read x-y CSV-style files
    """
    if ':' in data_path:
        run(['scp', '-rC', data_path, '/tmp/plotting.csv'])
        data_path = '/tmp/plotting.csv'

    x = []
    y = []
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

    return np.array(x), np.array(y)

def csvWriter(filename, x, y):
    with open(filename, 'w') as f:
        writer = csv.writer(f)
        writer.writerows(zip(x, y))

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

    return np.array(centers), np.array(areas)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('files', nargs='*', help='Radial velocity profile csv files. Flowrate profile in reverse mode.')
    ap.add_argument('-p', '--porosity-profile', required=True, help='Radial porosity profile csv file')
    ap.add_argument('-R', default=5.01e-4, type=float, help='Column radius')
    ap.add_argument('--nrad', default=5, type=int, help='Radial discretization')
    ap.add_argument('--shelltype', default='EQUIDISTANT', help='Radial discretization type')
    ap.add_argument('--reverse', action='store_true', default=False, help='Take flowrate profile as input and calculate velocity profile')
    args = ap.parse_args()
    print(vars(args))

    centers, areas = getCentersAndAreas(args.R, args.nrad, args.shelltype)

    r0, eps = readfile(args.porosity_profile)

    if args.reverse: 
        for f in args.files:
            r,flow = readfile(f)

            assert np.allclose(r, centers, 1e-03)
            assert np.allclose(r, r0, 1e-03)

            u = flow / (areas * eps)
             

            csvWriter(f'velocities_{f}_new.csv',r,u)
            print(f"Avg velocity = {sum(u)/len(u)}")
    else: 
        for f in args.files:
            r,u = readfile(f)

            assert np.allclose(r, centers, 1e-03)
            assert np.allclose(r, r0, 1e-03)

            # The velocity calculated is averaged over an infinitesimal volume, not calculated at a cross section. And the porosity is inherent to the geometry. 
            flowrate = areas * u * eps
            csvWriter(f'flowrates_{f}_new.csv',r,flowrate)
            print(sum(flowrate))

if __name__ == "__main__":
    main()
