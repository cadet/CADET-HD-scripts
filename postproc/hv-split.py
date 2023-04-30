#!/usr/bin/env python3

"""hv-split.py
calculate the holdup volumes for a split chromatogram (from XNS -> 2DGRM)

THIS SCRIPT IS NOT CORRECT. 
============================
Conceptually, that is.
It assumes that for the 2D GRM , theres is no mass transfer across the compartments/shells. i.e., 0 radial dispersion 
"""

import sympy as sp
from sympy.parsing.sympy_parser import parse_expr
import json
import math
import argparse
from addict import Dict
import subprocess
from ruamel.yaml import YAML
from pathlib import Path

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
    cn = [ (1 - elem / cin) for elem in c]
    holdup_num = trapz(cn, x=t) * math.pi * R**2 * u
    return holdup_num

def num_holdup_vol_flowrate(t, c, flowrate, cin):
    cn = [ (1 - elem / cin) for elem in c]
    holdup_num = trapz(cn, x=t) * flowrate
    return holdup_num

def ana_holdup_vol(v_i, v_b, eps, qinf):
    holdup_real = (v_i)  + (eps * v_b)  +  ( (1-eps) * v_b * qinf )
    return holdup_real

def ana_nonbind_holdup_vol(v_i, v_b, eps):
    holdup_real = (v_i)  + (eps * v_b)
    return holdup_real

def ana_solid_holdup_vol(v_i):
    holdup_real = (v_i)
    return holdup_real

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

    return centers, areas

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('chromatograms', nargs='*', help='Split chromatograms')
    ap.add_argument('-pp','--porosity-profile', help='csv style file with <centers>,<porosities>')
    ap.add_argument('-vp','--velocity-profile', help='csv style file with <centers>,<velocities>')
    args = Dict(vars(ap.parse_args()))

    nrad = len(args.chromatograms)

    _,epsc = readChromatogram(args.porosity_profile)

    ## NOTE: INLET velocity profile?
    _,u = readChromatogram(args.velocity_profile)

    # u = [ 2.09e-4 ] * nrad

    qmax   =  4.88
    ka     = 1.144
    kd     = 2.0e-3

    R = 5.01e-4
    h = 0.016
    epsb = 0.75
    cin = 7.14e-3

    qinf = qmax * ka / (ka*cin + kd)

    centers, areas = getCentersAndAreas(R, nrad, radial_disc_type='EQUIDISTANT')

    for chrom,area,por,vel in zip(args.chromatograms,areas,epsc,u): 
        t,c = readChromatogram(chrom)

        vol_total = area * h
        vol_int = por * vol_total
        vol_beads = vol_total - vol_int

        vol_holdup_ana = ana_holdup_vol(vol_int, vol_beads, epsb, qinf)
        vol_holdup_chrom = num_holdup_vol_flowrate(t,c, vel,cin)

        print(vol_holdup_ana, vol_holdup_chrom, vol_holdup_chrom/vol_holdup_ana, vol_holdup_ana/vol_holdup_chrom)



if __name__ == "__main__":
    main()
