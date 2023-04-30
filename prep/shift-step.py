#!/usr/bin/env python3

"""
Shift a step function to the right. Useful to calculate inlet profiles at a
distance from the actual inlet (like the bed beginning).

Doesn't account for actual dispersion, but introduces numerical dispersion, but
it's not that important in case we calculate the cout-cin integral as long as
the outlet is far away from the inlet
"""

import argparse
import numpy as np

def main():
    ap = argparse.ArgumentParser()

    ap.add_argument('-u', '--velocity', type=float, help='Inlet velocity')
    ap.add_argument('-c', '--inlet-concentration', type=float, help='Inlet concentration')
    ap.add_argument('-l', '--length', type=float, help='Length of void region')
    ap.add_argument('-T', '--end-time', type=float, help='End time of curve')
    ap.add_argument('-n', '--numpoints', type=int, default=500, help='Number of equidistant points in the curve')

    args = ap.parse_args()

    dt = args.length / args.velocity
    t = np.linspace(0, args.end_time, args.numpoints)
    y = np.heaviside(t-dt, 1) * args.inlet_concentration
    # y = np.heaviside(t-dt, 1)

    # np.savetxt('shifted.csv', np.stack([t,y],axis=1), delimiter=',')

if __name__ == "__main__":
    main()
