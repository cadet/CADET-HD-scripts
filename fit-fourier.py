#!/usr/bin/env python3

"""
Given a discrete x,y data file, compute the best fit fourier series.
Plot and write coefficients to file

I use this to try to get some information out of the radial velocity profile convergence with different meshes.
"""

import numpy as np
from subprocess import run
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
import argparse
from scipy.interpolate import make_interp_spline, BSpline

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

    return np.array(x), np.array(y), xticks

def series(An, r):
    ## Amplitude phase form as per wikipedia
    sum = np.zeros_like(r)
    sum += An[0] / 2
    period = r[-1]
    for n, an in enumerate(An[1:]): 
        sum += an*np.cos(2*np.pi*(n+1)*r/period)
        # sum += an*np.sin(2*np.pi*(n+1)*r/period)
    return sum


def residual(An, r, signal):
    return signal - series(An, r)

def csvWriter(filename, x, y):
    import csv
    with open(filename, 'w') as f:
        writer = csv.writer(f)
        writer.writerows(zip(x, y))

def driver(fname:str, Nu=0, smoothen_order=0, smoothen_npoints=250):
    x,y,_ = readfile(fname)

    if smoothen_order: 
        xsmooth = np.linspace(min(x), max(x), smoothen_npoints) 
        x_new = np.array(x)
        y_new = np.array(y)
        ordering = x_new.argsort()
        x_new = x_new[ordering]
        y_new = y_new[ordering]
        spl = make_interp_spline(x_new, y_new, k=smoothen_order)  # type: BSpline
        ysmooth = spl(xsmooth)
        x = xsmooth
        y = ysmooth


    if Nu == 0: 
        Nu = len(y)

    # Least-squares -------------------
    # fitting Fourier series
    An = np.arange(Nu)

    res = least_squares(residual, An, args=(x, y), method='trf')
    An = res.x

    print(f'{fname}, {max(An)}')

    with plt.style.context('science'):
        fig, ax_dict = plt.subplot_mosaic([['plot', 'coeffs']], constrained_layout=True, figsize=(7,2))

        ax1 = ax_dict['plot']
        ax1.plot(x/np.max(x), y/np.max(y), label="signal")
        ax1.plot(x/np.max(x), series(An, x) / np.max(y), label="series", ls='dashed')
        # ax1.set_xlabel('x')
        # ax1.set_ylabel('y')
        ax1.set_title(r'Signal \& Fourier Series')
        ax1.legend()

        ax2 = ax_dict['coeffs']
        # ax2.bar(range(Nu), An, width=0.5, label='coeffs')
        ax2.plot(range(Nu), An, label='coeffs')
        # ax2.set_xlabel('Coeff Index')
        # ax2.set_ylabel('Value')
        ax2.set_title(f'Coeffs, N={Nu}')

        # csvWriter(f"{fname}_coeffs.csv", range(Nu), An)

        # plt.savefig(f'{fname}.pdf')
        plt.show()


def main(): 
    ap = argparse.ArgumentParser()
    ap.add_argument('FILES', nargs='*', help='files to process. ideally, smooth')
    ap.add_argument('-n', type=int, default=0, help='Nu value = number of fourier series terms. If not given, assumed equal to number of points in data')
    ap.add_argument("--smoothen", nargs='?', type=int, const=3, help="(Cubic) spline interpolation through provided data. --smoothen <order>, with default order=3")
    ap.add_argument("--smoothen-npoints", default=250, type=int, help="Number of points to use for --smoothen interpolation")
    args = ap.parse_args()

    # TODO: Parallelize
    for fname in args.FILES: 
        driver(fname, Nu=args.n, smoothen_order=args.smoothen, smoothen_npoints=args.smoothen_npoints)

if __name__ == "__main__": 
    main()
