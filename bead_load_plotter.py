#!/bin/env python3

"""
@file: bead_load_plotter.py
@desc: plot chromatograms for individual beads. Use it after `vis.py -b`

vis.py can be invoked in parallel batches of, say, n timesteps each.
Then just cat the timestep data together and run this script on the three output files:
    - bead_loading.inf: containing nts, nbeads and ncv
    - bead_loading.xyzr: containing bead xyzr info
    - bead_loading.dat: containing the integrated concentration values in (nts,nbeads,ncv) dimensions

where
    nts: number of timesteps
    nbeads: number of beads
    ncv: number of colorvars (scalars)

Typically, ncv = 1 since we only need 'q' values, and we take care of this in the vis.py step: `vis.py -b -c scalar_1`

All 3 of the above files need to be in little endian binary formats (int, double, double).

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


ap = argparse.ArgumentParser()
# ap.add_argument("files", nargs='*', help="files to plot")
ap.add_argument("-f", "--filebasename", default="bead_loading",
        help="basename of files with .inf, .xyzr, and .dat to read data from. Default is bead_loading")
ap.add_argument("-t", "--title", required=False, default="Particle loading",
        help="title")
ap.add_argument("-x", "--xlabel", required=False, default="Time",
        help="xlabel")
ap.add_argument("-y", "--ylabel", required=False, default="Normalized concentration",
        help="ylabel")
ap.add_argument("-n", "--normalize", required=False, action='store_true',
        help="normalize y data to the last value")
ap.add_argument("-o", "--output", required=False,
        help="output file")
ap.add_argument("-s", "--sort", required=False, choices=['z', 'xy', 'r'],
        help="color lines by z, xy or r=sqrt(x^2+y^2)")
ap.add_argument("-xl", "--xlims", required=False,nargs=2, type=float,
        help="x axis limits")
ap.add_argument("-yl", "--ylims", required=False,nargs=2, type=float,
        help="y axis limits")
ap.add_argument("--timesteps", help="Timesteps file. ASCII, newline separated.")
args = vars(ap.parse_args())

info = bin_to_arr(args['filebasename'] +  '.inf', "=i")
nts = info[0]
nbeads = info[1]
ncv = info[2]
xyzr = np.array(bin_to_arr(args['filebasename'] + '.xyzr', "=d")).reshape((nbeads,4))
data = np.array(bin_to_arr(args['filebasename'] + '.dat', "=d")).reshape((nts,nbeads,ncv))

x = []
if args['timesteps']:
    with open(args['timesteps'], 'r') as inputfile:
        x = [float(time) for time in inputfile]
else:
    x = list(range(nts))


if args['sort'] == 'z':
    ordering = xyzr[:,2].argsort()
    data = data[:,ordering,:]
    xyzr = xyzr[ordering]
elif args['sort'] == 'r':
    ordering = xyzr[:,3].argsort()
    data = data[:,ordering,:]
    xyzr = xyzr[ordering]
elif args['sort'] == 'xy':
    x2y2 = np.array([ v[0]**2 + v[1]**2 for v in xyzr ])
    ordering = x2y2.argsort()
    data = data[:,ordering,:]
    xyzr = xyzr[ordering]

## TODO: np.graadient to see if the slope is zero yet

with plt.style.context(['science']):
    fig, ax = plt.subplots()
    xs = []
    ys = []
    lines = []
    count =0
    color = None

    for ibead,color in zip(range(nbeads), cm.rainbow(np.linspace(0,1,nbeads))):     #type: ignore
        if args['normalize']:
            y = normalize(data[:,ibead,0])
        else:
            y = data[:,ibead,0]
        ax.plot(x, y, c=color)
    ax.set(title=args['title'])
    ax.set(xlabel=args['xlabel'])
    ax.set(ylabel=args['ylabel'])
    ax.autoscale(tight=True)
    if args['xlims']:
        plt.xlim(args['xlims'])
    if args['ylims']:
        plt.ylim(args['ylims'])
    if args['output']:
        fig.savefig(args['output'], dpi=300)
    plt.show()
