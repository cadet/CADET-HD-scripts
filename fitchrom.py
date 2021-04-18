#!/usr/bin/env python3

"""
Fits monodisperse cadet with xns chromatogram.

@deps: cadet, cadet-python
"""

from scipy.interpolate import interp1d
from scipy.optimize import minimize
import subprocess
import argparse
from cadet import Cadet

cadetpath = "/home/jayghoshter/local/bin/cadet-cli"
Cadet.cadet_path = cadetpath

count = 0

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

def sse(x0, y0, x, y):
    f = interp1d(x0, y0)
    y0new = f(x)

    sse_value = sum([(n1 - n2)**2 for n1, n2 in zip(y, y0new)])
    return sse_value

def loadh5(col_dispersion, filename, ref_curve_filename):

    sim = Cadet()
    sim.filename = filename
    sim.load()

    sim.root.input.model.unit_002.col_dispersion = col_dispersion

    sim.save()
    runout = sim.run()
    if runout.returncode != 0:
        raise RuntimeError
    sim.load()

    x = sim.root.output.solution.solution_times
    y = sim.root.output.solution.unit_003.solution_outlet_comp_000

    x0, y0 = readChromatogram(ref_curve_filename)

    # NOTE: Ensure that section times do not exceed available reference time
    sim.root.input.solver.sections.section_times = [min(x0), max(x0)]

    sse_value = sse(x0,y0, x,y)
    global count
    count = count + 1
    print(count, ':', col_dispersion, sse_value)

    return(sse_value)

def main():

    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs='*', help="files to plot")
    ap.add_argument("-r", "--reference", required=True, help="Reference chromatogram to fit everything to")
    ap.add_argument("-i", "--init", type=float, default=6.5e-7, help="Initial guess for the dispersion coefficient")
    args = vars(ap.parse_args())

    for fitfile in args['files']:
        out = minimize(
                loadh5,
                args['init'],
                args=(fitfile, args['reference']),
                options={
                    'gtol':1e-6,
                    'eps':1e-12,
                    'maxiter': 50
                    }
                )

main()
