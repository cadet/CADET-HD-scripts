#!/bin/python3

# NOTE: depends on the mplstyle files from https://github.com/garrettj403/SciencePlots
# Usage: ./plotChrom.py <file> <file>
# Output: plot.pdf and plot.jpg

import matplotlib.pyplot as plt
import sys
import argparse

def normalize(data):
    return [ x/data[-1] for x in data ]

def rescale(data, factor):
    return [x*factor for x in data]

def readChromatogram(data_path, delimiter):
    time= []
    conc= []
    with open(data_path, newline='') as csvfile:
        # data = list(csv.reader(csvfile))
        for line in csvfile:
            data_line = line.strip().split(delimiter)
            data_line = list(filter(None, data_line))
            if (data_line != []):
                time.append(float(data_line[0]))
                conc.append(float(data_line[1]))
    return time, conc

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
ap.add_argument("-n", "--normalize", required=False, action='store_true',
        help="normalize y data to the last value")
ap.add_argument("-nl", "--no-legend", required=False, action='store_true',
        help="don't show legend")
ap.add_argument("-o", "--output", required=False,
        help="output file")
ap.add_argument("-r", "--rescale", type=float, required=False,
        help="Rescale graph with multiplicative factor (1/holdup-ratio)")
ap.add_argument("-sr", "--save-rescaled", required=False, action='store_true',
        help="File to save rescaled data")
ap.add_argument("--csv", required=False, action='store_true',
        help="input is csv file instead of default space separated")
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

if args['csv']:
    delimiter = ','
else:
    delimiter = ' '

with plt.style.context(['science']):
    fig, ax = plt.subplots()
    xs = []
    ys = []
    lines = []
    count =0
    for filename,label,linestyle,marker in zip(args['files'], args['labels'], args['linestyles'], args['markers']):
        x, y = readChromatogram(filename, delimiter)
        if args['normalize']:
            y = normalize(y)
        if args['rescale']:
            x = rescale(x, args['rescale'])
        line = ax.plot(x, y, label=label, marker=marker, linestyle=linestyle)
        if args['fill']:
            ax.fill_between(x, y, 1, interpolate=True)
        xs.append(x)
        ys.append(y)
        lines.append(line)
        count+=1
    if not args['no_legend']:
        legend = ax.legend(loc='best', shadow=True)
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
            with open(filename+'-postscaled', 'w') as f:
                writer = csv.writer(f, delimiter=delimiter)
                writer.writerows(zip(x, y))

    if args['to_csv']:
        import csv
        for filename,x,y in zip(args['files'], xs, ys):
            with open(filename+'.csv', 'w') as f:
                writer = csv.writer(f, delimiter=',')
                writer.writerows(zip(x, y))
