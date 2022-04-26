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
import random

from matplotlib.ticker import ScalarFormatter,AutoMinorLocator
from matplotlib.ticker import FormatStrFormatter, StrMethodFormatter
from matplotlib import ticker


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
ap.add_argument("-f", "--filebasename", default="bead_loading", help="basename of files with .inf, .xyzr, and .dat to read data from. Default is bead_loading")

ap.add_argument("-t", "--title", default="Particle loading", help="title")
ap.add_argument("-x", "--xlabel", default="Time", help="xlabel")
ap.add_argument("-y", "--ylabel", default="$\\frac{c_s}{c_s^{max}}$", help="ylabel")
ap.add_argument("-xl", "--xlims", nargs=2, type=float, help="x axis limits")
ap.add_argument("-yl", "--ylims", nargs=2, type=float, help="y axis limits")
ap.add_argument("-lw", "--linewidth", default=1, type=float, help="linewidth")
ap.add_argument("--alpha", default=1, type=float, help="opacity")

ap.add_argument("-n", "--normalize", action='store_true', help="normalize y data to the last value")

ap.add_argument("-s", "--sort", choices=['z', 'xy', 'r'], help="color lines by z, xy or r=sqrt(x^2+y^2)")

ap.add_argument("-o", "--output", help="output file")

ap.add_argument("--timesteps", help="Timesteps file. ASCII, newline separated.")
ap.add_argument("--scatter", action='store_true', help="Timesteps file. ASCII, newline separated.")

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

if args['scatter']: 
    with plt.style.context(['science']):
        fig, ax = plt.subplots(figsize=(4,3))
        xs = []
        ys = []
        lines = []
        count =0
        color = None

        min_dy = 1e-5
        bt = 0.9

        for ibead,color in zip(range(nbeads), cm.rainbow(np.linspace(0,1,nbeads))):     #type: ignore
            if args['normalize']:
                y = normalize(data[:,ibead,0])
            else:
                y = data[:,ibead,0]

            ## Breakthrough percentage based calc
            y = np.array(y)
            ind = np.where(y > bt * max(y))[0]
            ind = ind[0]
            xt = x[ind]

            # ## slope based calculation
            # dy = np.gradient(y, x)
            # # ax.plot(x, dy)
            # ind = np.where(dy < min_dy)[0]
            # ## WARNING: Hack to do away with initial zeroes
            # ind = ind[ind>10][0]
            # xt = x[ind]

            if args['sort'] == 'z':
                yt = xyzr[ibead, 2]
            elif args['sort'] == 'r':
                yt = xyzr[ibead, 3]
            elif args['sort'] == 'xy':
                yt = np.sqrt(xyzr[ibead,0] ** 2  + xyzr[ibead,1] ** 2)
            else: 
                raise RuntimeError("Unknown sort method")

            xs.append(xt)
            ys.append(yt)

        plt.scatter(ys, xs, s=1)
        ax.set(title=args['title'])
        # ax.set(ylabel=f"Time for slope $<$ {min_dy}")
        ax.set(ylabel="$t_{s}^{90}$")
        ax.set(xlabel=args['sort'])

        ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

        ax.autoscale(tight=True)
        if args['xlims']:
            plt.xlim(args['xlims'])
        if args['ylims']:
            plt.ylim(args['ylims'])


        if args['output']:
            fig.savefig(args['output'], dpi=300)
        else: 
            plt.show()

else: 
    with plt.style.context(['science']):
        fig, ax = plt.subplots(figsize=(4,3))
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

            ax.plot(x, y, c=color, lw=args['linewidth'], alpha=args['alpha'], zorder=random.randrange(nbeads))
        ax.set(title=args['title'])
        ax.set(xlabel=args['xlabel'])
        ax.set(ylabel=args['ylabel'])
        ax.autoscale(tight=True)

        ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

        if args['xlims']:
            plt.xlim(args['xlims'])
        if args['ylims']:
            plt.ylim(args['ylims'])
        if args['output']:
            fig.savefig(args['output'], dpi=300)
        else: 
            plt.show()
