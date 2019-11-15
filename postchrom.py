#!/usr/bin/env python3

"""
@desc: Script to calculate holdup volumes in XNS Chromatography simulations.
@inputs: chromatogram, xns.flow.in, xns.mass.in, mesh.log
@author: Jayghosh S. Rao

"""
# TODO: handle numbers in dicts properly
# TODO: Allow use of mixdinfo: for use with older sims for checking
# TODO: consistent meshvolume factors/magnitude
# TODO: Handle rngdexp****
# TODO: writeout mesh data

import numpy as np
import sys
import os
import json

EXTRACT_KEYS = [
        'Mesh Scaling Factor',
        'Cylinder Radius',
        'Cylinder Volume',
        'Real Int Volume',
        'Real Bead Volume',
        'Column Length',
        'rCyl',
        'Mesh.ScalingFactor',
        ' Mesh volume (physical -1 | dimension 3)',
        ' Mesh volume (physical 5 | dimension 3)',
        ' Mesh volume (physical 6 | dimension 3)' ]
PHYS_MAP = {
        "physical 5": "Mesh Int Volume",
        "physical 6": "Mesh Bead Volume",
        "physical -1": "Mesh Total Volume"
}


def read_chromatogram(data_path):
    time= []
    conc= []
    with open(data_path, newline='') as csvfile:
        # data = list(csv.reader(csvfile))
        for line in csvfile:
            data_line = line.strip().split(' ')
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

def num_holdup_vol(t, c, R, u):
    cn = [ (1 - elem / c[-1]) for elem in c]
    holdup_num = np.trapz(cn, x=t) * np.pi * R**2 * u
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

def read_mesh_log(data_path):
    mesh_data = {}
    with open(data_path, 'r') as mesh_log_file:
        for line in mesh_log_file:
            data_line = line.strip().split(':')
            data_line = list(filter(None, data_line))
            if (data_line != []) and (data_line[0][0] != '#'):
                ## To extract volume data (key: value)
                if data_line[0] in EXTRACT_KEYS:
                    mesh_data.update({data_line[0]: data_line[-1]})
                ## To extract info: <key> : <value>
                if (len(data_line)>1) and (data_line[1] in EXTRACT_KEYS):
                    key = [ PHYS_MAP[item] for item in PHYS_MAP if item in data_line[1] ]
                    mesh_data.update({key[0]: data_line[2]})

                ##to extract settings: Mesh.ScalingFactor
                ##delete when not needed
                data_line = line.strip().split(' ')
                data_line = list(filter(None, data_line))
                if data_line[0] in EXTRACT_KEYS:
                    mesh_data.update({data_line[0]: data_line[1]})

    # print(json.dumps(mesh_data, indent=4, sort_keys=True))

    return mesh_data

def scale_mesh_volumes(mesh_data):
    if ("Mesh Scaling Factor" in mesh_data):
        msf = float(mesh_data['Mesh Scaling Factor'])
    elif ("Mesh.ScalingFactor" in mesh_data):
        msf = float(mesh_data['Mesh.ScalingFactor'])
    else:
        print("Error finding mesh scaling factor!")
        sys.exit(-1)

    for key in PHYS_MAP:
        mkey = PHYS_MAP[key]
        value = str( float(mesh_data[mkey]) * (msf ** 3) )
        mesh_data.update({mkey:value})


def main():
    fsimfolder = sys.argv[1]
    msimfolder = sys.argv[2]
    simtype = sys.argv[3]


    for file in os.listdir():
        if file.endswith(".log"):
            mesh_log_path = (os.path.abspath(file))

    t, c          = read_chromatogram(os.path.join('MASS', msimfolder, 'chromatogram'))
    xns_flow_data = read_xns_in(os.path.join('FLOW', fsimfolder, 'xns.in'))
    xns_mass_data = read_xns_in(os.path.join('MASS', msimfolder, 'xns.in'))
    mesh_data     = read_mesh_log(mesh_log_path)
    scale_mesh_volumes(mesh_data)

    print(json.dumps(mesh_data, indent=4))

    try:
        eps  = float(xns_mass_data['clcepsilon'][0])
        qmax = float(xns_mass_data['mclqmax'][0])
        ka   = float(xns_mass_data['mclka'][0])
        kd   = float(xns_mass_data['mclkd'][0])
        cin  = float(xns_mass_data['rngdexp'][2])
    except:
        pass

    try:
        qinf = qmax * ka / (ka*cin + kd)
    except:
        pass

    u    = float(xns_flow_data['rngdexp'][2])

    # print(json.dumps(xns_flow_data, indent=4))
    # print(json.dumps(xns_mass_data, indent=4))


    R        = float(mesh_data['Cylinder Radius'])
    v_i_real = float(mesh_data['Real Int Volume'])
    v_b_real = float(mesh_data['Real Bead Volume'])
    v_i_mesh = float(mesh_data['Mesh Int Volume'])
    v_b_mesh = float(mesh_data['Mesh Bead Volume'])

    if (simtype == 'b'):
        vol_real_holdup = ana_holdup_vol(v_i_real, v_b_real, eps, qinf)
        vol_mesh_holdup = ana_holdup_vol(v_i_mesh, v_b_mesh, eps, qinf)
    elif (simtype == 'n'):
        vol_real_holdup = ana_nonbind_holdup_vol(v_i_real, v_b_real, eps)
        vol_mesh_holdup = ana_nonbind_holdup_vol(v_i_mesh, v_b_mesh, eps)
    elif (simtype == 's'):
        vol_real_holdup = ana_solid_holdup_vol(v_i_real)
        vol_mesh_holdup = ana_solid_holdup_vol(v_i_mesh)
    else:
        print("Bad simtype expression: use 's | n | b' ")

    vol_nume_holdup = num_holdup_vol(t, c, R, u)

    print("real Vh = ", vol_real_holdup)
    print("mesh Vh = ", vol_mesh_holdup)
    print("nume Vh = ", vol_nume_holdup)

    print("")
    n_m = vol_nume_holdup/vol_mesh_holdup
    n_r = vol_nume_holdup/vol_real_holdup
    print("nume/mesh Vh = ", round(n_m, 2)," = ", round((n_m-1)*100,2) , "%" )
    print("nume/real Vh = ", round(n_r, 2)," = ", round((n_r-1)*100,2) , "%" )

    print("")

if __name__ == "__main__":
    main()
