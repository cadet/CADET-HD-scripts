#!/usr/bin/env python3

"""
Apply a pointwise linear regression-based interpolation to predict the next curve in a series.

Input: Files corresponding to individual curves, extradimensional score associated with each curve

For ex. If I have data from a mesh convergence study, and want to predict the the next curve corresponding to a mesh size, the inputs are the files for the existing curves, and the scores would be the mesh sizes.

"""

import csv
import os
import numpy as np
from subprocess import run

from matplotlib import pyplot as plt
from sklearn.linear_model import LinearRegression
import argparse

def fit_line(x, y): 
    """
    Apply a linear regression to a given x,y data
    Return the R2 score, slope and intercept
    """
    X = np.array(x).reshape(-1,1)
    Y = np.array(y)

    model = LinearRegression()
    model.fit(X, Y)

    score = model.score(X,Y)

    return score, model.coef_[0], model.intercept_

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

    return x, y

def csvWriter(filename, x, y):
    with open(filename, 'w') as f:
        writer = csv.writer(f)
        writer.writerows(zip(x, y))

def main(): 
    ap = argparse.ArgumentParser()
    ap.add_argument('FILES', nargs='*', help='Files with 2-column csv data')
    ap.add_argument('--scores', nargs='*', type=float, required=True, help='scores corresponding to the input FILES')
    ap.add_argument('--predict', type=float, required=True, help='The score for which the curve will be predicted')
    args = ap.parse_args()

    ys = []
    xs = []
    final_y = []

    with plt.style.context('science'): 
        fig, ax = plt.subplots()

        for f in args.FILES: 
            x,y = readfile(f)
            ys.append(y)
            xs.append(x)
            ax.plot(x,y, lw=2)

        ys = np.array(ys)
        yst = np.transpose(ys) # (npoint, nmesh)

        n_points = len(ys[0])

        for i in range(n_points): 
            r2, m,c = fit_line(args.scores, yst[i])
            if r2 < 0.90: 
                print(f"Bad-ish fit for {i}: R2 = {r2}")
            final_y.append(m * args.predict + c)

        # WARNING: Assumes that x is the same
        ax.plot(xs[0], final_y, lw=2, c='black')
        plt.show()

        csvWriter(f'predicted-{args.predict}.csv', xs[0], final_y)

if __name__ == "__main__": 
    main()



