#!/usr/bin/env python3

"""
Apply a pointwise linear regression-based interpolation to predict the next curve in a series.

Input: Files corresponding to individual curves, extradimensional score associated with each curve

For ex. If I have data from a mesh convergence study, and want to predict the the next curve corresponding to a mesh size, the inputs are the files for the existing curves, and the scores would be the mesh sizes.

"""

import csv
import os
import numpy as np

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

def read_file(fname:str): 
    """
    Read a 2-column csv file
    """
    x = []
    y = []
    with open(fname, 'r') as fp:
        for line in csv.reader(fp): 
            x.append(float(line[0]))
            y.append(float(line[1]))

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
            x,y = read_file(f)
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



