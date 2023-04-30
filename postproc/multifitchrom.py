#!/usr/bin/env python3

"""
Based on fitchrom.py

Program to fit 2D GRM bulk data to XNS bulk data (bulk = interstitial)
"""

from scipy.interpolate import interp1d
from scipy.optimize import minimize
import argparse
from cadet import Cadet

cadetpath = "/home/jayghoshter/local/bin/cadet-cli"
Cadet.cadet_path = cadetpath

def sse2(x0, y0, x, y):
    """
    Sum of squares of errors
    x0, y0 are reference values
    x, y are to be measured against x0, y0
    """
    # Interpolate reference data (x0,y0) onto new x grid
    f = interp1d(x0, y0)
    y0new = f(x)

    sse_value = sum([(n1 - n2)**2 for n1, n2 in zip(y, y0new)])
    return sse_value

def sse(y0, y):
    sse_value = sum([(n1 - n2)**2 for n1, n2 in zip(y, y0)])
    return sse_value

def objective(params, h5file, target_bin_file):
    """ Objective function """

    ## Load cadet h5
    sim = Cadet()
    sim.filename = h5file
    sim.load()

    ## Modify parameters
    ## TODO:
    sim.root.input.model.unit_002.col_dispersion = params
    sim.root.input.model.unit_002.rad_dispersion = params

    # NOTE: Ensure that section times do not exceed available reference time
    sim.root.input.solver.sections.section_times = [min(x0), max(x0)]

    ## Save and run cadet
    sim.save()
    runout = sim.run()
    if runout.returncode != 0:
        print(runout)
        raise RuntimeError
    sim.load()

    y = sim.root.output.solution.unit_002.solution_bulk_comp_000

    y0 = load_reference_binary(target_bin_file)

    sse_value = sse(y0, y)
    print(sse_value)
