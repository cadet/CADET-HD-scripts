#!/usr/bin/env python3

"""
Script to evaluate the holdup volume per radial shell


"""

import math
import argparse
import numpy as np
import csv

eps  = 0.75
qmax = 4.88
ka   = 1.144
kd   = 2.0e-3
cin  = 7.14e-3

qinf = qmax * ka / (ka*cin + kd)

def csvWriter(filename, x, y):
    with open(filename, 'w') as f:
        writer = csv.writer(f)
        writer.writerows(zip(x, y))

def readfile(data_path, columns=[0,1], header=False):
    """ Read x-y CSV-style files
    """
    if ':' in data_path:
        run(['scp', '-rC', data_path, '/tmp/plotting.csv'])
        data_path = '/tmp/plotting.csv'

    x = []
    y = []
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
    return x, y

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
    
def ana_holdup_vol(v_i, v_b, eps, qinf):
    holdup_real = (v_i)  + (eps * v_b)  +  ( (1-eps) * v_b * qinf )
    return holdup_real

def num_holdup_vol(cin, cout, t, flowrate, eps):
    cn = [ (ein - eout)/ cin[-1] for ein,eout in zip(cin, cout) ]
    holdup_num = trapz(cn, x=t) * flowrate 
    return holdup_num

def trapz(y, *args, **kwargs):
     x = kwargs.get('x', range(len(y)))
     sum = 0.0
     for i in range(len(x) - 1):
         sum = sum + 0.5 * (y[i] + y[i+1]) * (x[i+1] - x[i])
     return sum

def main(): 
    ap = argparse.ArgumentParser()
    ap.add_argument('-p', '--porosity', required=True, help='radial porosity profile')
    ap.add_argument('-u', '--velocity', help='radial velocity profile')
    ap.add_argument('-f', '--flowrate', help='radial velocity profile')
    ap.add_argument('-ic', '--inlet-chromatograms', nargs='*', help='chromatograms to calculate the numerical holdup volumes')
    ap.add_argument('-oc', '--outlet-chromatograms', nargs='*', help='chromatograms to calculate the numerical holdup volumes')
    ap.add_argument('-R', '--column-radius', type=float, required=True, help='column radius')
    ap.add_argument('-L', '--column-length', type=float, required=True, help='column length')
    ap.add_argument('-st', '--shelltype', default= 'EQUIDISTANT', help='column radius')
    args = ap.parse_args()

    pcenters, porosity_profile = readfile(args.porosity)
    nrad = len(porosity_profile)
    centers, areas = getCentersAndAreas(args.column_radius, nrad, args.shelltype)

    if args.velocity: 
        vcenters, velocity_profile = readfile(args.velocity)
        assert np.allclose(pcenters, vcenters, rtol=1e-3)
        flowrates = [ u*A*e for u,A,e in zip(velocity_profile, areas, porosity_profile)]
    elif args.flowrate: 
        fcenters, flowrates = readfile(args.flowrate)
        assert np.allclose(pcenters, fcenters, rtol=1e-3)

    # assert len(porosity_profile) == len(velocity_profile)
    # np.allclose(pcenters, vcenters, rtol=1e-6)

    assert np.allclose(centers,  pcenters, rtol=1e-3)

    volumes = [ x * args.column_length for x in areas ]

    int_volumes = [ x * p for x,p in zip(volumes,porosity_profile) ]
    bead_volumes = [ x-y for x,y in zip(volumes, int_volumes) ]

    hv_ana = []
    hv_num = []
    ts = []
    cs = []

    for vi, vb in zip(int_volumes, bead_volumes): 
        hv_ana.append ( ana_holdup_vol(vi, vb, eps, qinf) )

    for chrom_in, chrom_out, flow, porosity in zip(args.inlet_chromatograms, args.outlet_chromatograms, flowrates, porosity_profile):
        tin,cin = readfile(chrom_in)
        tout,cout = readfile(chrom_out)

        ## Assume time series are equal
        assert tin == tout

        ts.append(tout)
        cs.append(cout)

        # WARNING: This REQUIRES us to either use (c_in-c_out) or (c_bulk + c_par) 
        # to calculate the correct holdup volume
        hv_num.append(num_holdup_vol(cin, cout, tin, flow, porosity))

    np.savetxt('hv_ana_rad.csv', np.stack([centers, hv_ana], axis=1), delimiter=',')
    np.savetxt('hv_num_rad.csv', np.stack([centers, hv_num], axis=1), delimiter=',')
    np.savetxt('hv_ratio_rad.csv', np.stack([centers, [x/y for x,y in zip(hv_num,hv_ana)]], axis=1), delimiter=',')

    ## Correct the given chromatograms using the calculated HV Ratios
    for fname,t,c,hva,hvn in zip(args.outlet_chromatograms,ts,cs,hv_ana,hv_num): 
        t = [ x * hva/hvn for x in t ]
        csvWriter(f'{fname}_corrected.csv', t, c)


if __name__ == "__main__":
    main()

