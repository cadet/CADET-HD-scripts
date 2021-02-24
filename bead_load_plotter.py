#!/bin/env python3

"""
@file: bead_load_plotter.py
@desc: plot chromatograms for individual beads. Use it after `vis.py -b`

vis.py can be invoked in parallel batches of, say, n timesteps each.
Then just cat the timestep data together and run this script on the three output files:
    - bead_loading.inf: containing nts, nbeads and ncv
    - bead_loading.xyzr: containing bead xyzr info
    - bead_loading.dat: containing the integrated concentration values in (nts,nbeads,ncv) dimensions

"""

from matplotlib import pyplot as plt
import matplotlib.cm as cm
import pickle
import numpy as np
import struct
import argparse

def bin_to_arr(filename, f):
    with(open(filename, 'rb')) as input:
        myiter = struct.iter_unpack(f, input.read())

        arr = []
        for i in myiter:
            arr.append(i[0])

        return arr

def unpickler(filename):
    with open(filename, "rb") as f:
        item = pickle.load(f)
        print("Unpickled", filename)
    return item

def normalize(data):
    return [ x/data[-1] for x in data ]

info = bin_to_arr('bead_loading.inf', "=i")
nts = info[0]
nbeads = info[1]
ncv = info[2]
xyzr = np.array(bin_to_arr('bead_loading.xyzr', "=d")).reshape((nbeads,4))
data = np.array(bin_to_arr('full.dat', "=d")).reshape((nts,nbeads,ncv))


ap = argparse.ArgumentParser()
ap.add_argument("files", nargs='*', help="files to plot")
ap.add_argument("-t", "--title", required=False, default="Chromatogram",
        help="title")
ap.add_argument("-x", "--xlabel", required=False, default="Time",
        help="xlabel")
ap.add_argument("-y", "--ylabel", required=False, default="Concentration",
        help="ylabel")
ap.add_argument("-n", "--normalize", required=False, action='store_true',
        help="normalize y data to the last value")
ap.add_argument("-o", "--output", required=False,
        help="output file")
ap.add_argument("-s", "--sort", required=False, choices=['z', 'xy', 'r'],
        help="output file")
args = vars(ap.parse_args())

if args['sort'] == 'z':
    xyzr = xyzr[xyzr[:,2].argsort()]
    data = data[:,xyzr[:,2].argsort(),:]
elif args['sort'] == 'r':
    xyzr = xyzr[xyzr[:,3].argsort()]
    data = data[:,xyzr[:,3].argsort(),:]
elif args['sort'] == 'xy':
    x2y2 = np.array([ v[0]**2 + v[1]**2 for v in xyzr ])
    xyzr = xyzr[x2y2.argsort()]
    data = data[:,x2y2.argsort(),:]

with plt.style.context(['science']):
    fig, ax = plt.subplots()
    xs = []
    ys = []
    lines = []
    count =0
    color = None

    for ibead,color in zip(range(nbeads), cm.rainbow(np.linspace(0,1,nbeads))):
        if args['normalize']:
            y = normalize(data[:,ibead,0])
        else:
            y = data[:,ibead,0]
        ax.plot(y, c=color)
    ax.set(title=args['title'])
    ax.set(xlabel=args['xlabel'])
    ax.set(ylabel=args['ylabel'])
    ax.autoscale(tight=True)
    if args['output']:
        fig.savefig(args['output'])
    plt.show()
