#!/bin/python3

# TODO: Output table in the right format
# TODO: Automatic handling of ticks, labels
# TODO: easy handling of time/speedup/efficiency outputs

import matplotlib.pyplot as plt
import csv
import sys


def readfile(data_path):
    x = []
    y = []
    with open(data_path, newline='') as csvfile:
        # data = list(csv.reader(csvfile))
        for line in csvfile:
            data_line = line.strip().split(' ')
            data_line = list(filter(None, data_line))
            if (data_line != []):
                x.append(float(data_line[0]))
                y.append(float(data_line[1]))
    return x, y


with plt.style.context(['science']):
    fig, ax = plt.subplots()
    xs = []
    ys = []
    speedups=[]
    lines = []
    count =0
    for arg in sys.argv[1:]:
        x, y = readfile(arg)
        speedup = [(item/y[0]) for item in y]
        # ideal = [(item/x[0]) for item in x]
        # efficiency = [s/i for s,i in zip(speedup,ideal)]
        efficiency = [(y[0]/item)*100 for item in y]

        # print("Cores Time Speedup Ideal Efficiency")
        # for i in range(len(x)):
        #     print(x[i], y[i], speedup[i], ideal[i], efficiency[i])

        # print(x)
        # print(speedup)
        # print(ideal)
        # print(efficiency)

        line = ax.plot(x, speedup, label=arg, marker='s')
        # line = ax.plot(x, efficiency, label=arg, marker='s')
        lines.append(line)
        count+=1
    # lineref = ax.plot(x, ideal, label='ideal', linestyle='--')
    # lines.append(lineref)
    legend = ax.legend(loc='best', shadow=True)
    ax.set_xscale('log')
    # ax.set_yscale('log')
    ax.set(title='Speedup Plot')
    ax.set(xlabel='Cores')
    ax.set(ylabel='Speedup')
    plt.minorticks_off()
    ax.autoscale(tight=True)
    plt.xticks([120, 240, 480, 960, 1920, 3840, 4100], ['120', '240', '480', '960', '1920', '3840', ''])
    # plt.yticks([1, 2, 4, 8, 16, 32], ['1', '2', '4', '8', '16', '32'])
    fig.savefig(sys.argv[1]+'.pdf')
    fig.savefig(sys.argv[1]+'.jpg', dpi=300)
