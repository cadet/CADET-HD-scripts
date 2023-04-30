#!/usr/bin/env python3

"""
Calculate the ideal total holdup volume and divide the difference equally among the radial sections
"""

import argparse
import numpy as np

eps  = 0.75
qmax = 4.88
ka   = 1.144
kd   = 2.0e-3
cin  = 7.14e-3
qinf = qmax * ka / (ka*cin + kd)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('files', nargs='*', help='Files with radial output chromatograms to transform')
    ap.add_argument('-f', '--flowrates', help='flowrates')
    ap.add_argument('-cin', '--inlet_concentration', type=float, help='concentration value')

    ap.add_argument('-p', '--porosity', type=float, help='Column porosity value')
    ap.add_argument('-L', '--length', type=float, help='Column Length')
    ap.add_argument('-R', '--radius', type=float, help='Column Radius')


    args = ap.parse_args()

    flowrates = np.loadtxt(args.flowrates, delimiter=',').T[1]

    num_hvs = []
    arrs = []

    for f, flow in zip(args.files, flowrates): 
        arr = np.loadtxt(f, delimiter=',').T
        arrs.append(arr)
        num_hvs.append( np.trapz( 1.0 - arr[1]/args.inlet_concentration  , x=arr[0]) * flow )

    print(sum(num_hvs))

    total_volume = np.pi * args.radius**2 * args.length
    int_volume = args.porosity * total_volume
    bed_volume = total_volume - int_volume

    ana_hv = ana_holdup_vol(int_volume, bed_volume, eps, qinf)
    print(ana_hv)

    delta = sum(num_hvs) - ana_hv
    nrad = len(flowrates)
    print(f'{nrad = }')

    new_hvs = []

    for f,arr,flow in zip(args.files,arrs,flowrates): 
        # ratio = (hvn-delta/nrad)/hvn
        ratio = ana_hv/sum(num_hvs)
        t_new = arr[0] * ratio
        new_hvs.append( np.trapz( 1.0 - arr[1]/args.inlet_concentration  , x=t_new) * flow )
        np.savetxt(f'{f}_corrected.csv', np.stack([t_new, arr[1]], axis=1), delimiter=',')
    print(sum(new_hvs))

def ana_holdup_vol(v_i, v_b, eps, qinf):
    holdup_real = (v_i)  + (eps * v_b)  +  ( (1-eps) * v_b * qinf )
    return holdup_real

if __name__ == "__main__":
    main()
