#!/bin/python3

# depends on the mplstyle files from https://github.com/garrettj403/SciencePlots
# Usage: ./plotChrom.py <file> <file>
# Output: plot.pdf and plot.jpg

# TODO: argparse support for xlabel, ylabel and title
# TODO: argparse support for plt.show
# TODO: argparse support for delimiter

import matplotlib.pyplot as plt
# import csv
import sys

def readChromatogram(data_path):
    time= []
    conc= []
    with open(data_path, newline='') as csvfile:
        # data = list(csv.reader(csvfile))
        for line in csvfile:
            data_line = line.strip().split(' ')
            data_line = list(filter(None, data_line))
            if (data_line != []):
                time.append(float(data_line[0]))
                conc.append(float(data_line[1]))
    return time, conc


with plt.style.context(['science']):
    fig, ax = plt.subplots()
    xs = []
    ys = []
    lines = []
    count =0
    for arg in sys.argv[1:]:
        x, y = readChromatogram(arg)
        line = ax.plot(x, y, label=arg)
        xs.append(x)
        ys.append(y)
        lines.append(line)
        count+=1
    # legend = ax.legend(loc='best', shadow=True)
    ax.set(title='Chromatogram')
    ax.set(xlabel='Time')
    ax.set(ylabel='Concentration')
    ax.autoscale(tight=True)
    fig.savefig('plot.pdf')
    fig.savefig('plot.jpg', dpi=300)
