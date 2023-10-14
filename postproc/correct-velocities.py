#!/usr/bin/env python3

"""
Given a velocity profile and a target total flowrate, correct the velocity profile by shifting it up/down such that the total flowrate is matched. 

XNS interstitial velocity -> Corrected.
"""

import argparse
import numpy as np
import math
from scipy.optimize import minimize

def getCentersAndAreas(col_radius:float, nrad:int, radial_disc_type='EQUIDISTANT'): 
    nRegions = nrad
    nShells = nRegions + 1 #Including r = 0
    rShells = []

    ## NOTE: move to current unit?
    # R = params.col_radius
    R = col_radius

    if radial_disc_type == 'EQUIVOLUME':
        for n in range(nShells):
            rShells.append(R * math.sqrt(n/nRegions))
    elif radial_disc_type == 'EQUIDISTANT':
        for n in range(nShells):
            rShells.append(R * (n/nRegions))

    centers = [ (rShells[i+1] + rShells[i])/2 for i in range(nRegions) ]
    areas = [ math.pi *  (rShells[i+1]**2 - rShells[i]**2) for i in range(nRegions) ]

    return np.array(centers), np.array(areas)

def get_total_flowrate(ueps:float, areas, u, p): 
    return sum(areas * p * (u + ueps))

def get_total_flowrate_error(ueps:float, areas, u, p, f0): 
    return abs(get_total_flowrate(ueps, areas, u, p) - f0)/f0


def main(): 
    ap = argparse.ArgumentParser()
    ap.add_argument('velocity_profiles', nargs='*', help='files with velocity profiles')
    ap.add_argument('-tf', '--target-flowrate', type=float, help='Target total flowrate for the final velocity profile')
    ap.add_argument('-p', '--porosity-profiles', nargs='*', help='corresponding porosity profiles' )
    ap.add_argument('-R', '--column-radius', type=float, help='column radius')
    ap.add_argument('-st', '--shelltype', choices=['EQUIDISTANT', 'EQUIVOLUME'], default='EQUIDISTANT', help='type of shell discretization')
    ap.add_argument('--nrad', type=int, help='number of radial shells')
    
    args = ap.parse_args()

    centers, areas = getCentersAndAreas(args.column_radius, args.nrad, args.shelltype)

    for ufile, pfile in zip(args.velocity_profiles, args.porosity_profiles): 
        print(ufile, pfile)

        u = np.loadtxt(ufile, delimiter=',').T
        p = np.loadtxt(pfile, delimiter=',').T

        assert np.allclose(u[0], p[0], rtol=1e-3)

        out = minimize(
                get_total_flowrate_error,
                0.0,
                args=(areas, u[1], p[1], args.target_flowrate),
                options={
                    'gtol':1e-12,
                    'eps':1e-12,
                    'maxiter': 1000
                    }
                )

        print(out.x)
        u_new = np.array(u)
        u_new[0] = u[0]
        u_new[1] = u[1] + out.x

        np.savetxt(ufile + '_corrected.csv', u_new.T, delimiter=',' )

if __name__ == "__main__":
    main()
