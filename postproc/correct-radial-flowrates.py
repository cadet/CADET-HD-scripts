#!/usr/bin/env python

"""
Pushes flowrate profiles up or down to achieve a target total flowrate. Thus I can calculate the XNS flowrate, and use this to remove the flowrate error while preserving the 2DGRM flowrate profile. Ideally, I use XNS velocity -> forced flowrate -> Corrected forced flowrates.
"""

import argparse
import numpy as np

def main(): 
    ap = argparse.ArgumentParser()
    ap.add_argument('files', nargs='*', help='list of flowrate profiles to correct')
    ap.add_argument('-t', '--target', type=float, help='target flowrate to achieve')
    ap.add_argument('-o', '--output', help='output filename')
    args = ap.parse_args()

    for f in args.files: 
        arr = np.loadtxt(f, delimiter=',').T
        print(arr[1])

        delta = args.target - sum(arr[1])
        nrad = arr.shape[1]

        arr[1] = arr[1] + delta/nrad
        print(arr[1])

        ofname = args.output or f'{f}_corrected.csv'

        np.savetxt(ofname, arr.T, delimiter=',')
    

if __name__ == "__main__": 
    main()
