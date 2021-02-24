#!/bin/python3

# NOTE: depends on the mplstyle files from https://github.com/garrettj403/SciencePlots
# Usage: ./plotChrom.py <file> <file>
# Output: plot.pdf and plot.jpg

import matplotlib.pyplot as plt
import sys
import argparse

def normalize(data, refValue):
    return [ x/refValue for x in data ]

def rescale(data, factor):
    return [x*factor for x in data]

def rescaleVariant(x, y, factor):
    import numpy as np
    factors = np.linspace(1,factor,len(x))
    # factors = [ 1 + (factor-1) * j/max(y) for j in y]
    return [ i*j for i,j in zip(x,factors) ]
    # return [ i*factor*(j/max(y)) for i,j in zip(x,y) ]

# def readChromatogram(data_path, delimiter):
def readChromatogram(data_path):
    time= []
    conc= []
    delimiter = ' '
    with open(data_path, newline='') as csvfile:
        if ',' in csvfile.readline():
            delimiter = ','
    with open(data_path, newline='') as csvfile:
        # data = list(csv.reader(csvfile))
        for line in csvfile:
            data_line = line.strip().split(delimiter)
            data_line = list(filter(None, data_line))
            if (data_line != []):
                time.append(float(data_line[0]))
                conc.append(float(data_line[1]))
    return time, conc

def trapz(y, *args, **kwargs):
     x = kwargs.get('x', range(len(y)))
     sum = 0.0
     for i in range(len(x) - 1):
         sum = sum + 0.5 * (y[i] + y[i+1]) * (x[i+1] - x[i])
     return sum

def num_holdup_vol(t, c, R, u):
    import math
    cn = [ (1 - elem / c[-1]) for elem in c]
    holdup_num = trapz(cn, x=t) * math.pi * R**2 * u
    return holdup_num

ap = argparse.ArgumentParser()
ap.add_argument("files", nargs='*', help="files to plot")
ap.add_argument("-t", "--title", required=False, default="Chromatogram",
        help="title")
ap.add_argument("-x", "--xlabel", required=False, default="Time",
        help="xlabel")
ap.add_argument("-y", "--ylabel", required=False, default="Concentration",
        help="ylabel")
ap.add_argument("-l", "--labels", required=False, nargs='*',
        help="legend labels")
ap.add_argument("-ls", "--linestyles", required=False, nargs='*',
        help="linestyles = solid dashed ...")
ap.add_argument("-m", "--markers", required=False, nargs='*',
        help="markers = s, o, ...")
ap.add_argument("-xl", "--xlims", required=False,nargs=2, type=float,
        help="x axis limits")
ap.add_argument("-yl", "--ylims", required=False,nargs=2, type=float,
        help="y axis limits")
ap.add_argument("-f", "--fill", required=False, action='store_true',
        help="fill area under curve")
ap.add_argument("-n", "--normalize", required=False,
        help="normalize y data to the reference value")
ap.add_argument("-nl", "--no-legend", required=False, action='store_true',
        help="don't show legend")
ap.add_argument("--legend", required=False, nargs=3, default=['upper center', '0.5', '-0.2'], type=str,
        help="Legend settings: --legend <location> <bbox_to_anchor>")
ap.add_argument("-o", "--output", required=False,
        help="output file")
ap.add_argument("-r", "--rescale", type=float, required=False,
        help="Rescale graph with multiplicative factor (1/holdup-ratio)")
ap.add_argument("-rv", "--rescale-variant", type=float, required=False,
        help="Rescale graph based on y value with multiplicative factor (1/holdup-ratio)")
ap.add_argument("-sr", "--save-rescaled", required=False,
        help="File to save rescaled data")
ap.add_argument("-hv", "--holdup-volume", required=False, nargs = 2, type=float,
        help="File to save rescaled data")
# ap.add_argument("--csv", required=False, action='store_true',
#         help="input is csv file instead of default space separated")
ap.add_argument("--to-csv", required=False, action='store_true',
        help="save as csv")
args = vars(ap.parse_args())

if not args['labels']:
    args['labels'] = args['files']

if not args['linestyles']:
    args['linestyles'] = [ 'solid' ] * len(args['files'])
    # args['linestyles'] = []
elif len(args['linestyles']) == 1:
    args['linestyles'] = args['linestyles'] * len(args['files'])

if not args['markers']:
    args['markers'] = [ None ] * len(args['files'])
    # args['markers'] = []
elif len(args['markers']) == 1:
    args['markers'] = args['markers'] * len(args['files'])

# if args['csv']:
#     delimiter = ','
# else:
#     delimiter = ' '

with plt.style.context(['science']):
    fig, ax = plt.subplots()
    xs = []
    ys = []
    lines = []
    count =0
    for filename,label,linestyle,marker in zip(args['files'], args['labels'], args['linestyles'], args['markers']):
        # x, y = readChromatogram(filename, delimiter)
        x, y = readChromatogram(filename)
        if args['normalize']:
            y = normalize(y, args['normalize'])
        if args['rescale']:
            x = rescale(x, args['rescale'])
        if args['rescale_variant']:
            x = rescaleVariant(x, y, args['rescale_variant'])
        line = ax.plot(x, y, label=label, marker=marker, linestyle=linestyle)
        if args['fill']:
            ax.fill_between(x, y, 1, interpolate=True)
        xs.append(x)
        ys.append(y)
        lines.append(line)
        count+=1
    if not args['no_legend']:
        # legend = ax.legend(loc='best', shadow=True)
        legend = ax.legend(loc=args['legend'][0], bbox_to_anchor=(float(args['legend'][1]),float(args['legend'][2])), shadow=True)
    ax.set(title=args['title'])
    ax.set(xlabel=args['xlabel'])
    ax.set(ylabel=args['ylabel'])
    ax.autoscale(tight=True)
    if args['xlims']:
        plt.xlim(args['xlims'])
    if args['ylims']:
        plt.ylim(args['ylims'])
    if args['output']:
        fig.savefig(args['output'])
    plt.show()


    if args['save_rescaled']:
        import csv
        for filename,x,y in zip(args['files'], xs, ys):
            with open(filename+'-' + args['save_rescaled'], 'w') as f:
                writer = csv.writer(f, delimiter=',')
                writer.writerows(zip(x, y))

    if args['normalize']:
        import csv
        for filename,x,y in zip(args['files'], xs, ys):
            with open(filename+'-normalized', 'w') as f:
                writer = csv.writer(f, delimiter=',')
                writer.writerows(zip(x, y))


    if args['to_csv']:
        import csv
        for filename,x,y in zip(args['files'], xs, ys):
            with open(filename+'.csv', 'w') as f:
                writer = csv.writer(f, delimiter=',')
                writer.writerows(zip(x, y))

    if args['holdup_volume']:
        for filename,x,y in zip(args['files'], xs, ys):
            holdup = num_holdup_vol(x, y, args['holdup_volume'][0], args['holdup_volume'][1])
            print("{filename}: {holdup}".format(filename=filename, holdup=holdup))


