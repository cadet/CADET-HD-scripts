#!/usr/bin/env python3

""" postchrom2.py

Successor to postchrom.py

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
import gmsh

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

def read_xns_in(data_path):
    xns_dict = {}
    with open(data_path, 'r') as xns_in_file:
        for line in xns_in_file:
            data_line = line.strip().split(' ')
            data_line = list(filter(None, data_line))
            if (data_line != []) and (data_line[0][0] != '#'):
                xns_dict.update({data_line[0] : data_line[1:]})
    return(xns_dict)

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

def ana_holdup_vol(v_i, v_b, eps, qinf):
    holdup_real = (v_i)  + (eps * v_b)  +  ( (1-eps) * v_b * qinf )
    return holdup_real

def ana_nonbind_holdup_vol(v_i, v_b, eps):
    holdup_real = (v_i)  + (eps * v_b)
    return holdup_real

def ana_solid_holdup_vol(v_i):
    holdup_real = (v_i)
    return holdup_real

def process_rngdexp(f,R):
    x, y, r, theta = sp.symbols('x y r theta')
    # f = 2.0*2.09e-4*( 1.0 - (x*x + y*y)/( 5.00999444e-4 * 5.00999444e-4 ) )
    f = parse_expr(f)
    f_r = f.subs({x: r*sp.cos(theta), y: r*sp.sin(theta)}).simplify()
    f_r_2pr = (f_r*2*r*(math.pi)).simplify()
    int_val=sp.integrate(f_r_2pr, (r, 0, R))
    # act_int_val = math.pi * R * R * 2.09e-4
    avg = int_val/(math.pi * R * R)
    print("Average rngdexp = ", avg)
    return float(avg)

def gmsh_calculate_volumes(filename, nthreads:int = 1): 
    print('GMSH API:', gmsh.GMSH_API_VERSION)

    gmsh_version = 'unknown'
    try:
        gmsh_version = subprocess.check_output(["gmsh", "--version"], stderr=subprocess.STDOUT).strip().decode('utf8')
    except subprocess.CalledProcessError:
        pass

    print('GMSH Version:', gmsh_version)

    gmsh.initialize()

    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("General.NumThreads", nthreads)

    gmsh.merge(filename)

    print("[Column]")
    gmsh.plugin.setNumber("MeshVolume", "PhysicalGroup", -1)
    gmsh.plugin.setNumber("MeshVolume", "Dimension", 3)
    t = gmsh.plugin.run("MeshVolume")

    _, _, data = gmsh.view.getListData(t)
    volume_column = data[0][-1]

    print("[Interstitial]")
    gmsh.plugin.setNumber("MeshVolume", "PhysicalGroup", 5)
    gmsh.plugin.setNumber("MeshVolume", "Dimension", 3)
    t = gmsh.plugin.run("MeshVolume")

    _, _, data = gmsh.view.getListData(t)
    volume_interstitial = data[0][-1]

    print("[Packed Bed]")
    gmsh.plugin.setNumber("MeshVolume", "PhysicalGroup", 6)
    gmsh.plugin.setNumber("MeshVolume", "Dimension", 3)
    t = gmsh.plugin.run("MeshVolume")

    _, _, data = gmsh.view.getListData(t)
    volume_packedBed = data[0][-1]

    return {'column': volume_column, 'interstitial': volume_interstitial, 'bed': volume_packedBed}

def main():

    ap = argparse.ArgumentParser()
    ap.add_argument('-nt', '--nthreads', default=1, type=int)
    ap.add_argument('-xf','--xns-flow-in')
    ap.add_argument('-xm','--xns-mass-in')
    ap.add_argument('-c','--chromatogram')
    ap.add_argument('-t', '--type', default='binding', choices = ['binding', 'nonbinding', 'solid'])

    ap.add_argument('-m','--mesh')
    ap.add_argument('-vi','--volume-interstitial', type=float)
    ap.add_argument('-vb','--volume-bed', type=float)

    ap.add_argument('-mc','--mesh-config')
    ap.add_argument('-R', '--column-radius', type=float)

    args = Dict(vars(ap.parse_args()))

    simtype = args.type

    t, c          = readChromatogram(args.chromatogram)

    flow_in = Path(args.xns_flow_in) 
    mass_in = Path(args.xns_mass_in)

    if flow_in.is_dir(): 
        flow_in = flow_in / 'xns.in'

    if mass_in.is_dir(): 
        mass_in = mass_in / 'xns.in'

    xns_flow_data = read_xns_in(flow_in)
    xns_mass_data = read_xns_in(mass_in)

    if args.volume_interstitial and args.volume_bed: 
        mesh_volumes = {
                'interstitial': args.volume_interstitial,
                'bed': args.volume_bed
                }
    else: 
        mesh_volumes = gmsh_calculate_volumes(args.mesh, args.nthreads)


    try:
        eps  = float(xns_mass_data['clcepsilon'][0])
        qmax = float(xns_mass_data['mclqmax'][0])
        ka   = float(xns_mass_data['mclka'][0])
        kd   = float(xns_mass_data['mclkd'][0])
        cin  = float(xns_mass_data['rngdexp'][2])
    except:
        raise RuntimeError("ERROR: Couldn't read xns mass input file!")

    try:
        qinf = qmax * ka / (ka*cin + kd)
    except:
        raise RuntimeError("ERROR: Couldn't calculate qinf!")

    # print(json.dumps(xns_flow_data, indent=4))
    # print(json.dumps(xns_mass_data, indent=4))


    if args.column_radius: 
        R = args.column_radius
    else: 
        # Read mesh_config
        mesh_config = Dict(YAML(typ='safe').load(Path(args.mesh_config)))

        msf     = mesh_config.gmsh['Mesh.ScalingFactor']
        R       = mesh_config.container.size[6] * msf

        if mesh_config.container.shape != 'cylinder': 
            raise NotImplementedError

    # h        = mesh_config.container.size[5] * msf

    # v_full_ideal = math.pi * R**2 * h
    # v_b_ideal = 
    # v_i_ideal = 

    v_i_mesh = mesh_volumes['interstitial']
    v_b_mesh = mesh_volumes['bed']

    u = xns_flow_data['rngdexp'][2:]
    u = ''.join(u)
    u_avg = process_rngdexp(u, R)

    if (simtype == 'binding'):
        # vol_ideal_holdup = ana_holdup_vol(v_i_ideal, v_b_ideal, eps, qinf)
        vol_mesh_holdup = ana_holdup_vol(v_i_mesh, v_b_mesh, eps, qinf)
    elif (simtype == 'nonbinding'):
        # vol_ideal_holdup = ana_nonbind_holdup_vol(v_i_ideal, v_b_ideal, eps)
        vol_mesh_holdup = ana_nonbind_holdup_vol(v_i_mesh, v_b_mesh, eps)
    elif (simtype == 'solid'):
        # vol_ideal_holdup = ana_solid_holdup_vol(v_i_ideal)
        vol_mesh_holdup = ana_solid_holdup_vol(v_i_mesh)
    else:
        raise RuntimeError("Bad type")

    vol_nume_holdup = num_holdup_vol(t, c, R, u_avg, cin)

    n_m = vol_nume_holdup/vol_mesh_holdup
    # n_r = vol_nume_holdup/vol_ideal_holdup
    print("num/mesh Vh = ", round(n_m, 2)," = ", round((n_m-1)*100,2) , "%" )
    # print("num/ideal Vh = ", round(n_r, 2)," = ", round((n_r-1)*100,2) , "%" )
   
    holdup_volumes = { 'numerical': vol_nume_holdup, 'mesh': vol_mesh_holdup }
    holdup_volume_ratios = {'num_mesh': n_m}
    holdup_volume_errors = {'num_mesh': (n_m-1)*100 }

    result = Dict()
    result.volumes = mesh_volumes
    result.holdup_volumes = holdup_volumes
    result.holdup_volume_ratios = holdup_volume_ratios
    result.holdup_volume_errors = holdup_volume_errors

    with open('post.json', 'w') as fp:
        json.dump(result, fp, indent=4)
    
if __name__ == "__main__":
    main()
