#!/usr/bin/env python3

"""
USAGE: ./fluxes.py -R 5.01e-4 -u 2.09e-4 -c 7.14e-3 VC-mono-0.10-rngout.cg -o test.pdf
"""

import argparse
import math
import matplotlib.pyplot as plt 

def csvWriter(filename, x, y):
    import csv
    with open(filename, 'w') as f:
        writer = csv.writer(f)
        writer.writerows(zip(x, y))

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

def num_holdup_vol(t, c, R, u, cin):
    import math
    cn = [ (1 - elem / cin) for elem in c]
    holdup_num = trapz(cn, x=t) * math.pi * R**2 * u
    return holdup_num

def inlet_flux_total(time, R, u, c_in):
    return c_in * time * R**2 * math.pi * u

def outlet_flux_total(time, R, u, conc_array):
    return trapz(conc_array, x=time) * math.pi * R**2 * u

ap = argparse.ArgumentParser()
ap.add_argument("files", nargs='*', help="files to plot")
ap.add_argument("-R", "--column-radius", type=float, help="Column radius to calculate area")
ap.add_argument("-u", "--velocity", type=float, help="Flow velocity")
ap.add_argument("-c", "--concentration", type=float, help="Inlet concentration")
ap.add_argument("-o", "--output", help="output plot file")
args = vars(ap.parse_args())

R = args['column_radius']
u = args['velocity']
c_in = args['concentration']

with plt.style.context(['science']):
    fig, ax = plt.subplots()
    for filename in args['files']:
        t, c = readChromatogram(filename)

        influx = [ inlet_flux_total(current_time, R, u, c_in) for current_time in t ]
        outflux = [ outlet_flux_total(t[0:i], R, u, c[0:i]) for i in range(len(c))]
        stored = [ x-y for x,y in zip(influx,outflux)]

        ax.plot(t, influx)
        ax.plot(t, outflux)

        if args['output']:
            fig.savefig(args['output'])
        else:
            plt.show()
        csvWriter(f'{filename}_influx.csv', t, influx)
        csvWriter(f'{filename}_outflux.csv', t, outflux)
        csvWriter(f'{filename}_stored.csv', t, stored)

