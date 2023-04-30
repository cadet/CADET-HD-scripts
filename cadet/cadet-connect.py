#!/usr/bin/env python3
"""
Create a CADET connections matrix given a set of flowrates/velocity
"""

import argparse 
import math 
import csv
from ruamel.yaml import YAML
from pathlib import Path
from subprocess import run
from rich import print

def readfile(data_path, columns=[0,1], header=False, xticksColumn=0):
    """ Read x-y CSV-style files
    """
    if ':' in data_path:
        run(['scp', '-rC', data_path, '/tmp/plotting.csv'])
        data_path = '/tmp/plotting.csv'

    x = []
    y = []
    xticks = []
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
                if xticksColumn is not None: 
                    xticks.append(float(data_line[xticksColumn]))

    return x, y, xticks

def csvWriter(filename, x, y):
    with open(filename, 'w') as f:
        writer = csv.writer(f)
        writer.writerows(zip(x, y))

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


def getFlowrates(nrad, velocity, radial_disc_type, col_radius):
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

    flowrates = [velocity * math.pi *  (rShells[i+1]**2 - rShells[i]**2) for i in range(nRegions)]
    return flowrates

def connectParallel(unit_port_left:list, unit_port_right:list, flowrates:list):
    """ Connect two units in parallel
    unit_port_left and unit_port_right are lists of tuples: [(0,0), (1,0), (2,0)] and [(3,0), (3,1), (3,2)]
    """
    connections = []
    assert(len(unit_port_left) == len(unit_port_right))

    for left, right, flowrate in zip(unit_port_left, unit_port_right, flowrates):
        line = [ left[0], right[0], left[1], right[1], -1, -1, flowrate]
        connections.extend(line)

    return connections

def connectSerial(unit_port_left:list, unit_port_right:list, flowrates:list):
    """ Connect two units serially 
    unit_port_left and unit_port_right are lists of tuples: [(3,0), (3,1), (3,2)] and [(4,0)]
    """
    connections = []

    pairs = [(x,y) for x in unit_port_left for y in unit_port_right]

    # print(unit_port_left)
    # print(unit_port_right)
    # print(pairs)
    # print(len(pairs))
    # print(len(flowrates))

    assert(len(pairs) == len(flowrates))

    for pair, flowrate in zip(pairs, flowrates):
        line = [pair[0][0], pair[1][0], pair[0][1], pair[1][1], -1, -1, flowrate]
        connections.extend(line)

    return connections


if __name__ == "__main__": 

    ap = argparse.ArgumentParser()
    ap.add_argument('-nr', '--nrad', type=int, default=1, help='radial discretization zones')
    ap.add_argument('-v','-u', '--velocity', type=float, help='velocity')
    ap.add_argument('-R', '--column-radius', type=float, help='column radius')
    ap.add_argument('-st', '--shelltype', choices = ['EQUIDISTANT', 'EQUIVOLUME'], default='EQUIDISTANT', help='Type of shell discretization')
    ap.add_argument('-c', '--connection', choices = ['1d', '1d_dpfr', '2d_serial_inlet', '2d_parallel_inlet', '2d_serial_dpfr' ], help='connection scheme')
    ap.add_argument('-f', '--flowrates', help='use flowrates from 2nd column of given csv file instead of using constant velocity from --velocity')
    args = ap.parse_args()

    print(vars(args))


    flowrates = []
    if args.flowrates: 
        print(f"[bold yellow]Using column 2 of file:{args.flowrates} as provided flowrates![/]")
        _,flowrates,_ = readfile(args.flowrates)
        assert(args.nrad == len(flowrates))

        if args.velocity and args.column_radius: 
            flowrates_constant_vel = getFlowrates(args.nrad, args.velocity, args.shelltype, args.column_radius)
            centers, areas = getCentersAndAreas(args.column_radius, args.nrad, args.shelltype)

            total_flowrate = sum(flowrates)
            print(f"Total flowrate = {total_flowrate}")

            expected_flowrate = args.velocity * math.pi * (args.column_radius)**2
            print(f"Expected total flowrate = {expected_flowrate}")

            ratio = total_flowrate / expected_flowrate
            print(f"Flowrate ratio = {ratio}")

            error = (1 - ratio) * 100
            print(f"Error = {error}%")

    else: 
        flowrates = getFlowrates(args.nrad, args.velocity, args.shelltype, args.column_radius) 

    connections = []

    if args.connection == '1d': 
        # CONNECTION: inlet -> 1d_column
        # flowrates = getFlowrates(1, args.velocity, args.shelltype, args.column_radius) 

        conn1 = connectSerial([(0,0)], [(1,0)], [sum(flowrates)])

        connections.extend(conn1)

    elif args.connection == '1d_dpfr': 
        # CONNECTION: inlet -> DPFR -> 1d_column -> DPFR
        # flowrates = getFlowrates(1, args.velocity, args.shelltype, args.column_radius) 

        conn1 = connectSerial([(0,0)], [(1,0)], flowrates)
        conn2 = connectSerial([(1,0)], [(2,0)], flowrates)
        conn3 = connectSerial([(2,0)], [(3,0)], flowrates)

        connections.extend(conn1)
        connections.extend(conn2)
        connections.extend(conn3)

    elif args.connection == '2d_parallel_inlet': 
        ## NOTE: Assumes unit_000 to unit_{nrad-1} are inlets, and unit_{nrad} is column
        # CONNECTION: (nrad)*inlets -> 2d_column -> outlet

        # flowrates = getFlowrates(args.nrad, args.velocity, args.shelltype, args.column_radius) 

        conn1 = connectParallel([(x,0) for x in range(args.nrad)], [(args.nrad, y) for y in range(args.nrad)], flowrates)
        conn2 = connectSerial([(args.nrad, y) for y in range(args.nrad)], [(args.nrad+1,0)], flowrates) 
        connections.extend(conn1)
        connections.extend(conn2)

    elif args.connection == '2d_serial_dpfr': 
        # CONNECTION: inlet -> DPFR -> 2d_column -> DPFR
        # flowrates = getFlowrates(1, args.velocity, args.shelltype, args.column_radius) 
        conn1 = connectSerial([(0,0)], [(1,0)], [sum(flowrates)])

        # flowrates = getFlowrates(args.nrad, args.velocity, args.shelltype, args.column_radius) 
        conn2 = connectSerial([(1,0)], [(2,y) for y in range(args.nrad)], flowrates)
        conn3 = connectSerial([(2,y) for y in range(args.nrad)], [(3,0)], flowrates)

        connections.extend(conn1)
        connections.extend(conn2)
        connections.extend(conn3)

    elif args.connection == '2d_serial_inlet': 
        # CONNECTION: inlet -> 2d_column -> outlet 

        conn1 = connectSerial([(0,0)], [(1,y) for y in range(args.nrad)], flowrates)
        conn2 = connectSerial([(1,y) for y in range(args.nrad)], [(2,0)], flowrates)

        connections.extend(conn1)
        connections.extend(conn2)

    YAML(typ='safe').dump({'connections': connections}, Path('connections.yaml'))
    print(f"{connections = }")
