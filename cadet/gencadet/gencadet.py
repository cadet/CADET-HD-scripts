#!/usr/bin/env python3

"""gencadet.py

A rewrite of xcad.py.
Hacked together from xcad.py, cadet-connect.py and yaml2hdf
Sometimes rewrites are quicker than refactors.

The idea is that I should be able to store all the information from a given XNS 3D simulation (reduced) into one file and be able to generate cadet simulations out of it.

So a mono xns simulation -> info 1d and 2d (using pack.py and paravision) --> yaml file
Then yaml file --> m1 and m2

EXPERIMENTAL / WORK IN PROGRESS
"""

# WARNING: Porosities have to be manually adjusted when changing between dpfr and non-dpfr modes. We don't have separate fields for this yet.
# WARNING: Only works for GRM and MULTI_COMPONENT_LANGMUIR currently

# TODO: Couple with pack.py? -> Porosities

from addict import Dict
import numpy as np
import argparse
from cadet import Cadet
from subprocess import run

from ruamel.yaml import YAML
from pathlib import Path
import csv
import math
from rich import print

cadetpath = run(['which', 'cadet-cli'], capture_output=True, text=True ).stdout.strip()
Cadet.cadet_path = cadetpath
ENCODING='ascii'
yaml=YAML(typ='safe')
# yaml.default_flow_style=False

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

def getFlowratesInterstitial(velocities:list, porosities:list, col_radius:float, nrad:int, radial_disc_type:str = 'EQUIDISTANT' ): 
    """
    Return flowrates corresponding to the given interstitial velocities.

    Interstitial velocities can be calculated from XNS. To enforce them in
    cadet, we need to convert them to a flowrate for the connections matrix.

    In 2DGRM the continuity equation is only respected per channel. So A*v is
    only conserved within a given channel. So the interstitial velocities
    calculated in CADET for a given inlet velocity will not match the
    interstitial velocities in XNS.

    In XNS, we see velocity directly proportional to porosity.
    In CADET, we see the opposite due to A*v_s = A*v_i*eps -> v_i = v_s/eps.

    So if we were to enforce the right interstitial velocities, we need to
    input them as extrapolated flowrates as we do here. 

    """
    assert len(velocities) == nrad
    assert len(porosities) == nrad
    _, areas = getCentersAndAreas(col_radius, nrad, radial_disc_type)
    flowrates = [ v * A * e for v,A,e in zip(velocities, areas, porosities) ] 
    return flowrates

def connectParallel(unit_port_left:list, unit_port_right:list, flowrates:list):
    """ Connect two units in parallel
    unit_port_left and unit_port_right are lists of tuples: [(0,0), (1,0), (2,0)] and [(3,0), (3,1), (3,2)]
    """
    connections = []
    assert(len(unit_port_left) == len(unit_port_right))
    assert(len(unit_port_left) == len(flowrates))

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


def get_inlet_unit(params):

    InletUnit                     = Dict()
    InletUnit.inlet_type          = 'PIECEWISE_CUBIC_POLY'
    InletUnit.unit_type           = 'INLET'
    InletUnit.ncomp               = params.ncomp
    InletUnit.sec_000.const_coeff = params.inlet_concentration

    InletUnit.ports = 1

    return InletUnit

def get_outlet_unit(params):

    OutletUnit                     = Dict()
    OutletUnit.unit_type           = 'OUTLET'
    OutletUnit.ncomp               = params.ncomp

    OutletUnit.ports = 1

    return OutletUnit

def get_dpfr_unit(params):
    DPFRUnit=Dict()

    DPFRUnit.unit_type = 'LUMPED_RATE_MODEL_WITHOUT_PORES'
    DPFRUnit.adsorption_model = 'NONE'

    DPFRUnit.ncomp = params.ncomp
    DPFRUnit.col_dispersion = 0
    DPFRUnit.col_length = params.void_length
    DPFRUnit.cross_section_area = np.pi * params.col_radius**2

    DPFRUnit.velocity = 1
    DPFRUnit.total_porosity = 1
    DPFRUnit.init_c = [0.0]
    DPFRUnit.init_q = [0.0]

    DPFRUnit.discretization.ncol = params.ncol_dpfr

    DPFRUnit.discretization.nbound = [0] * params.ncomp
    DPFRUnit.discretization.reconstruction = 'WENO'
    DPFRUnit.discretization.use_analytic_jacobian = 1
    DPFRUnit.discretization.weno.boundary_model = 0
    DPFRUnit.discretization.weno.weno_eps = 1e-10
    DPFRUnit.discretization.weno.weno_order = 3

    DPFRUnit.ports = 1

    return DPFRUnit

def get_lrm_unit(params):
    LRMUnit = Dict()

    LRMUnit.unit_type = 'LUMPED_RATE_MODEL_WITHOUT_PORES'

    LRMUnit.adsorption_model      = params.adsorption_model
    LRMUnit.adsorption.is_kinetic = params.is_kinetic

    ## Linear model
    LRMUnit.adsorption.lin_ka     = params.lin_ka
    LRMUnit.adsorption.lin_kd     = params.lin_kd

    LRMUnit.ncomp = params.ncomp
    LRMUnit.col_dispersion = params.col_dispersion
    LRMUnit.col_length = params.bed_length
    LRMUnit.cross_section_area = np.pi * params.col_radius**2

    LRMUnit.velocity = 1
    LRMUnit.total_porosity = params.col_porosity
    LRMUnit.init_c = params.init_c
    LRMUnit.init_q = params.init_q

    LRMUnit.discretization.ncol = params.ncol

    LRMUnit.discretization.nbound = [1] * params.ncomp
    LRMUnit.discretization.reconstruction = 'WENO'
    LRMUnit.discretization.use_analytic_jacobian = 1
    LRMUnit.discretization.weno.boundary_model = 0
    LRMUnit.discretization.weno.weno_eps = 1e-10
    LRMUnit.discretization.weno.weno_order = 3

    LRMUnit.ports = 1

    return LRMUnit

def get_grm_unit(params, psd=None):
    GRMUnit = Dict()

    GRMUnit.unit_type         = 'GENERAL_RATE_MODEL'
    GRMUnit.ncomp             = params.ncomp
    GRMUnit.col_dispersion    = params.col_dispersion
    GRMUnit.col_length        = params.bed_length
    GRMUnit.col_porosity      = params.col_porosity
    GRMUnit.film_diffusion    = params.film_diffusion
    GRMUnit.par_diffusion     = params.par_diffusion
    GRMUnit.par_porosity      = params.par_porosity
    GRMUnit.par_radius        = params.par_radius
    GRMUnit.par_surfdiffusion = params.par_surfdiffusion
    GRMUnit.init_c            = params.init_c
    GRMUnit.init_cp           = params.init_cp
    GRMUnit.init_q            = params.init_q

    GRMUnit.velocity           = 1
    GRMUnit.cross_section_area = np.pi * params.col_radius**2

    GRMUnit.adsorption_model      = params.adsorption_model
    GRMUnit.adsorption.is_kinetic = params.is_kinetic

    # ## Linear model
    # GRMUnit.adsorption.lin_ka     = params.lin_ka
    # GRMUnit.adsorption.lin_kd     = params.lin_kd

    ## MCL
    GRMUnit.adsorption.mcl_ka     = params.mcl_ka
    GRMUnit.adsorption.mcl_kd     = params.mcl_kd
    GRMUnit.adsorption.mcl_qmax   = params.mcl_qmax

    # ## SMA
    # GRMUnit.adsorption.sma_ka         = params.sma_ka
    # GRMUnit.adsorption.sma_kd         = params.sma_kd
    # GRMUnit.adsorption.sma_nu         = params.sma_nu
    # GRMUnit.adsorption.sma_sigma      = params.sma_sigma
    # GRMUnit.adsorption.sma_lambda     = params.sma_lambda
    # GRMUnit.adsorption.sma_refc0      = params.sma_refc0
    # GRMUnit.adsorption.sma_refq       = params.sma_refq

    GRMUnit.discretization.ncol                  = params.ncol
    GRMUnit.discretization.npar                  = params.npar

    GRMUnit.discretization.use_analytic_jacobian = 1
    GRMUnit.discretization.nbound                = [1] * params.ncomp

    GRMUnit.discretization.par_disc_type       = 'EQUIDISTANT_PAR'
    GRMUnit.discretization.schur_safety        = 1.0e-8
    GRMUnit.discretization.weno.boundary_model = 0
    GRMUnit.discretization.weno.weno_eps       = 1e-10
    GRMUnit.discretization.weno.weno_order     = 3
    GRMUnit.discretization.gs_type             = 1
    GRMUnit.discretization.max_krylov          = 0
    GRMUnit.discretization.max_restarts        = 10

    GRMUnit.ports = 1

    if psd == True:
        # par_radius = unpickler(os.path.expanduser('~/fullbedavgRads.pickle'))
        # par_type_volfrac = unpickler(os.path.expanduser('~/fullbedvolfrac.pickle'))
        npartype = len(params.par_radius_psd)
        GRMUnit.adsorption_model               = [ params.adsorption_model ]
        GRMUnit.adsorption_model_multiplex     = 1
        GRMUnit.discretization.par_disc_type   = ['EQUIDISTANT_PAR']
        GRMUnit.discretization.npar            = [ params.npar ] * npartype

        GRMUnit.par_radius                     = params.par_radius_psd
        GRMUnit.par_type_volfrac               = params.par_type_volfrac

    return GRMUnit

def get_2dgrm_unit(params, psd=None):
    # NOTE: We use col_porosity_2d and par_type_volfrac_2d because one json
    # file can hold all the information about a 3D column and we can simply run
    # xcad.py <file> --runtype
    # FIXME: However, this only allows switching between mono1d, mono2d Or poly1d and poly2d
    # Since the porosity profile for mono and poly are different, we cannot store them in the same file this way
    # Meaning: Cannot use one complete file for mono2d and poly2d. 
    # On second thought, the mono and poly columns are different columns altogether, they would need 2 separate files to specify them anyway

    GRM2DUnit = Dict()

    GRM2DUnit.unit_type                = 'GENERAL_RATE_MODEL_2D'
    GRM2DUnit.ncomp                    = params.ncomp
    GRM2DUnit.col_dispersion           = params.col_dispersion
    GRM2DUnit.col_dispersion_radial    = params.col_dispersion_radial
    GRM2DUnit.col_length               = params.bed_length
    # GRM2DUnit.col_porosity             = params.col_porosity
    GRM2DUnit.col_porosity             = params.col_porosity_2d
    GRM2DUnit.film_diffusion           = params.film_diffusion
    GRM2DUnit.par_diffusion            = params.par_diffusion
    GRM2DUnit.par_porosity             = params.par_porosity
    GRM2DUnit.par_radius               = params.par_radius
    GRM2DUnit.par_surfdiffusion        = params.par_surfdiffusion
    GRM2DUnit.init_c                   = params.init_c
    GRM2DUnit.init_cp                  = params.init_cp
    GRM2DUnit.init_q                   = params.init_q

    GRM2DUnit.velocity           = 1
    GRM2DUnit.cross_section_area = np.pi * params.col_radius**2

    GRM2DUnit.col_radius            = params.col_radius

    GRM2DUnit.adsorption_model                = [ params.adsorption_model ]
    GRM2DUnit.adsorption_model_multiplex      = 1
    GRM2DUnit.adsorption.is_kinetic           = params.is_kinetic

    GRM2DUnit.adsorption.mcl_ka               = params.mcl_ka
    GRM2DUnit.adsorption.mcl_kd               = params.mcl_kd
    GRM2DUnit.adsorption.mcl_qmax             = params.mcl_qmax

    GRM2DUnit.adsorption.lin_ka               = params.lin_ka
    GRM2DUnit.adsorption.lin_kd               = params.lin_kd

    nrad     = len(params.col_porosity_2d)
    GRM2DUnit.discretization.ncol                  = params.ncol
    GRM2DUnit.discretization.nrad                  = nrad
    GRM2DUnit.discretization.npar                  = params.npar
    GRM2DUnit.discretization.radial_disc_type      = params.radial_disc_type
    # GRM2DUnit.discretization.radial_disc_type      = 'EQUIVOLUME'
    GRM2DUnit.discretization.use_analytic_jacobian = 1
    GRM2DUnit.discretization.nbound                = [1] * params.ncomp
    GRM2DUnit.discretization.par_disc_type         = 'EQUIDISTANT_PAR'

    GRM2DUnit.discretization.schur_safety          = 1.0e-8
    GRM2DUnit.discretization.weno.boundary_model   = 0
    GRM2DUnit.discretization.weno.weno_eps         = 1e-10
    GRM2DUnit.discretization.weno.weno_order       = 3
    GRM2DUnit.discretization.gs_type               = 1
    GRM2DUnit.discretization.max_krylov            = 0
    GRM2DUnit.discretization.max_restarts          = 10


    GRM2DUnit.ports = nrad

    if psd:
        # par_radius = unpickler(os.path.expanduser('~/fullbedavgRads.pickle'))
        # par_type_volfrac = unpickler(os.path.expanduser('~/fullbedvolfrac.pickle'))
        npartype = len(params.par_radius_psd)
        GRM2DUnit.adsorption_model               = [ params.adsorption_model ]
        GRM2DUnit.adsorption_model_multiplex     = 1
        GRM2DUnit.discretization.par_disc_type   = ['EQUIDISTANT_PAR']
        GRM2DUnit.discretization.npar            = [ params.npar ] * npartype
        GRM2DUnit.discretization.nbound           = [1] * npartype

        GRM2DUnit.par_radius                     = params.par_radius_psd
        # GRM2DUnit.par_type_volfrac               = params.par_type_volfrac
        GRM2DUnit.par_type_volfrac               = params.par_type_volfrac_2d
        GRM2DUnit.par_type_volfrac_multiplex     = 1
        GRM2DUnit.par_porosity                   = [ params.par_porosity ] * npartype

    return GRM2DUnit

def get_ioreturn():
    IOReturn                                   = Dict()
    IOReturn.write_solution_times              = 1
    IOReturn.split_components_data             = 1
    IOReturn.split_ports_data                  = 1
    IOReturn.unit.write_sens_bulk              = 0
    IOReturn.unit.write_sens_flux              = 0
    IOReturn.unit.write_sens_inlet             = 0
    IOReturn.unit.write_sens_outlet            = 0
    IOReturn.unit.write_sens_particle          = 0
    IOReturn.unit.write_solution_bulk          = 1
    IOReturn.unit.write_solution_flux          = 1
    IOReturn.unit.write_solution_inlet         = 1
    IOReturn.unit.write_solution_outlet        = 1
    IOReturn.unit.write_solution_particle      = 1
    IOReturn.unit.write_solution_solid         = 1
    IOReturn.unit.write_sens_column            = 0
    IOReturn.unit.write_sens_column_inlet      = 0
    IOReturn.unit.write_sens_column_outlet     = 0
    IOReturn.unit.write_solution_column        = 0
    IOReturn.unit.write_solution_column_inlet  = 1
    IOReturn.unit.write_solution_column_outlet = 1
    return IOReturn


def get_flowrates_wrapper(params, ndim): 
    flowrates = None
    nrad     = len(params.col_porosity_2d) if ndim == 2 else 1

    if ndim == 1: 
        if params.velocity_interstitial: 
            print("Enforcing interstitial velocity!")
            flowrates = getFlowratesInterstitial([params.velocity_interstitial], [params.col_porosity], params.col_radius, nrad, params.radial_disc_type)
        elif params.flowrates: 
            print("Using provided flowrates!")
            flowrates = params.flowrates
        else: 
            print("Calculating flowrates based on constant inlet velocity!")
            flowrates = getFlowrates(nrad, params.velocity, params.radial_disc_type, params.col_radius) 
    elif ndim == 2: 
        if params.velocity_interstitial_2d: 
            print("Enforcing interstitial velocity!")
            flowrates = getFlowratesInterstitial(params.velocity_interstitial_2d, params.col_porosity_2d, params.col_radius, nrad, params.radial_disc_type)
        elif params.flowrates_2d: 
            print("Using provided flowrates!")
            flowrates = params.flowrates_2d
        else: 
            print("Calculating flowrates based on constant inlet velocity!")
            flowrates = getFlowrates(nrad, params.velocity, params.radial_disc_type, params.col_radius) 
    else: 
        raise RuntimeError("Bad ndim!")

    print(f"{sum(flowrates) = }")
    return flowrates


def create_simulation(mode, ndim, psd, params, outfile):
    model = Cadet()
    model.filename = outfile

    root                     = model.root

    root.input['return'] = get_ioreturn()

    nrad     = len(params.col_porosity_2d) if ndim == 2 else 1

    flowrates = get_flowrates_wrapper(params, ndim)

    connections = []

    if mode == '1d': 
        # CONNECTION: inlet -> 1d_column
        # flowrates = getFlowrates(1, params.velocity, params.shelltype, params.column_radius) 
        root.input.model.nunits  = 2
        root.input.model.unit_000 = get_inlet_unit(params)
        root.input.model.unit_001 = get_grm_unit(params, psd=psd)

        root.input['return'].unit_000 = get_ioreturn().unit
        root.input['return'].unit_001 = get_ioreturn().unit

        conn1 = connectSerial([(0,0)], [(1,0)], [sum(flowrates)])

        connections.extend(conn1)

    elif mode == '1d_dpfr': 
        # CONNECTION: inlet -> DPFR -> 1d_column -> DPFR
        # flowrates = getFlowrates(1, params.velocity, params.shelltype, params.column_radius) 

        root.input.model.nunits   = 4
        root.input.model.unit_000 = get_inlet_unit(params)
        root.input.model.unit_001 = get_dpfr_unit(params)
        root.input.model.unit_002 = get_grm_unit(params, psd)
        root.input.model.unit_003 = get_dpfr_unit(params)

        root.input['return'].unit_000 = get_ioreturn().unit
        root.input['return'].unit_001 = get_ioreturn().unit
        root.input['return'].unit_002 = get_ioreturn().unit
        root.input['return'].unit_003 = get_ioreturn().unit

        conn1 = connectSerial([(0,0)], [(1,0)], flowrates)
        conn2 = connectSerial([(1,0)], [(2,0)], flowrates)
        conn3 = connectSerial([(2,0)], [(3,0)], flowrates)

        connections.extend(conn1)
        connections.extend(conn2)
        connections.extend(conn3)

    elif mode == '2d_parallel_inlet': 
        ## NOTE: Assumes unit_000 to unit_{nrad-1} are inlets, and unit_{nrad} is column
        # CONNECTION: (nrad)*inlets -> 2d_column -> outlet

        # flowrates = getFlowrates(params.nrad, params.velocity, params.shelltype, params.column_radius) 

        nrad     = len(params.col_porosity_2d)
        root.input.model.nunits   = nrad + 2
        for index in range(nrad):
            root.input.model[f"unit_{index:03d}"] = get_inlet_unit(params)
            root.input['return'][f"unit_{index:03d}"] = get_ioreturn().unit

        root.input.model[f"unit_{nrad:03d}"] = get_2dgrm_unit(params, psd)
        root.input['return'][f"unit_{nrad:03d}"] = get_ioreturn().unit

        root.input.model[f"unit_{nrad+1:03d}"] = get_outlet_unit(params)
        root.input['return'][f"unit_{nrad+1:03d}"] = get_ioreturn().unit


        conn1 = connectParallel([(x,0) for x in range(nrad)], [(nrad, y) for y in range(nrad)], flowrates)
        conn2 = connectSerial([(nrad, y) for y in range(nrad)], [(nrad+1,0)], flowrates) 
        connections.extend(conn1)
        connections.extend(conn2)

    elif mode == '2d_serial_dpfr': 
        # CONNECTION: inlet -> DPFR -> 2d_column -> DPFR
        # flowrates = getFlowrates(1, params.velocity, params.shelltype, params.column_radius) 
        nrad     = len(params.col_porosity_2d)
        root.input.model.nunits   = 4
        root.input.model.unit_000 = get_inlet_unit(params)
        root.input.model.unit_001 = get_dpfr_unit(params)
        root.input.model.unit_002 = get_2dgrm_unit(params, psd)
        root.input.model.unit_003 = get_dpfr_unit(params)

        root.input['return'].unit_000 = get_ioreturn().unit
        root.input['return'].unit_001 = get_ioreturn().unit
        root.input['return'].unit_002 = get_ioreturn().unit
        root.input['return'].unit_003 = get_ioreturn().unit

        conn1 = connectSerial([(0,0)], [(1,0)], [sum(flowrates)])

        # flowrates = getFlowrates(params.nrad, params.velocity, params.shelltype, params.column_radius) 
        conn2 = connectSerial([(1,0)], [(2,y) for y in range(nrad)], flowrates)
        conn3 = connectSerial([(2,y) for y in range(nrad)], [(3,0)], flowrates)

        connections.extend(conn1)
        connections.extend(conn2)
        connections.extend(conn3)

    elif mode == '2d_serial_inlet': 
        # CONNECTION: inlet -> 2d_column -> outlet 
        nrad     = len(params.col_porosity_2d)
        root.input.model.nunits   = 3
        root.input.model.unit_000 = get_inlet_unit(params)
        root.input.model.unit_001 = get_2dgrm_unit(params, psd)
        root.input.model.unit_002 = get_outlet_unit(params)

        root.input['return'].unit_000 = get_ioreturn().unit
        root.input['return'].unit_001 = get_ioreturn().unit
        root.input['return'].unit_002 = get_ioreturn().unit
        root.input['return'].unit_003 = get_ioreturn().unit

        conn1 = connectSerial([(0,0)], [(1,y) for y in range(nrad)], flowrates)
        conn2 = connectSerial([(1,y) for y in range(nrad)], [(2,0)], flowrates)

        connections.extend(conn1)
        connections.extend(conn2)

    root.input.model.connections.nswitches                 = 1
    root.input.model.connections.switch_000.section        = 0
    root.input.model.connections.switch_000.connections    = connections
    root.input.model.connections.connections_include_ports = 1

    root.input.model.solver.gs_type       = 1
    root.input.model.solver.max_krylov    = 0
    root.input.model.solver.max_restarts  = 10
    root.input.model.solver.schur_safety  = 1e-8

    root.input.solver.sections.section_times         = params.section_times
    root.input.solver.user_solution_times            = np.linspace(0, params.section_times[-1], params.nts)
    root.input.solver.sections.nsec                  = len(params.section_times)-1
    root.input.solver.sections.section_continuity    = [0]* (len(params.section_times) - 2)
    root.input.solver.time_integrator.abstol         = 1e-10
    root.input.solver.time_integrator.algtol         = 1e-12
    root.input.solver.time_integrator.init_step_size = 1e-6
    root.input.solver.time_integrator.max_steps      = 1000000
    root.input.solver.time_integrator.reltol         = 1e-10
    root.input.solver.nthreads                       = 0

    return model

def simulate(mode, ndim, psd, params, outfile):
    sim = create_simulation(mode, ndim, psd, params, outfile)
    sim.save()
    save_yaml(reconstruct(sim.root), outfile.replace('.h5', '.yaml'))
    print(sim.run())
    sim.load()

    import csv
    for csvout in params.csvout:
        csvoutarr = csvout.split('/')
        x = sim.root.output.solution.solution_times
        y = sim.root.output.solution
        for item in csvoutarr:
            y = y[item]
        with open(outfile + '.' + csvout.replace('/', '_') +'.csv', 'w') as f:
            writer = csv.writer(f, delimiter=',')
            writer.writerows(zip(x, y))

    return sim

def read_simulation(infile):
    sim = Cadet()
    sim.filename = infile
    sim.load()
    return sim

## TODO: Might not be required
def reconstruct(dic:dict):
    """
    Convert from addict dict filled with numpy types to
    dict filled with python native types
    """
    # cleandict = {}
    dic = dict(dic)
    for key in dic.keys():
        if isinstance(dic[key], Dict):
            value = reconstruct(dic[key])
            # cleandict.update({key: value})
            dic.update({key: value})
        else:
            value = np2native(dic[key])
            # cleandict.update({key: value})
            dic.update({key: value})
    # return cleandict
    return dic

def update_nested(refdict:dict, newdict:dict):
    """
    Update a dictionary in a nested manner.
    (Don't just overwrite top level key/values as is default behavior)
    """
    for key, value in newdict.items():
        if key in refdict.keys():
            if isinstance(value, dict):
                update_nested(refdict[key], newdict[key])
            else:
                refdict[key] = value
        else:
            refdict.update({key: value})

def np2native(obj):
    """
    Convert from numpy types to python native types
    """
    if isinstance(obj, bytes):
        return obj.decode(ENCODING)
    if isinstance(obj, np.bytes_):
        return obj.tobytes().decode(ENCODING)
    elif isinstance(obj,np.ndarray):
        if any(isinstance(x, bytes) for x in obj):
            return [ x.decode(ENCODING) for x in obj ]
        elif any(isinstance(x, np.bytes_) for x in obj):
            return [ x.tobytes().decode(ENCODING) for x in obj ]
        else:
            return obj.tolist()
    elif isinstance(obj, np.generic):
        return obj.tolist()
    else:
        return obj

def load_file(fname, key=''):
    ext = Path(fname).suffix
    config = {}
    if ext == ".h5":
        config = load_h5(fname)
    elif ext == ".yaml" or ext == ".yml" :
        config = yaml.load(Path(fname))

    if key == 'input':
        config = delete_key(config, 'output')
    elif key == 'output':
        config = delete_key(config, 'input')

    return config

def save_file(config, fname, update=False):
    ext = Path(fname).suffix
    if ext == ".h5":
        out = save_h5(config, fname, update)
    elif ext == ".yaml" or ext == ".yml" :
        out = save_yaml(config, fname)
    return out # type: ignore

def save_yaml(config, fname):
    try:
        del(config['meta'])
    except KeyError:
        pass

    yaml.dump(config, Path(fname))

    return config

def load_h5(fname):
    sim = Cadet()
    sim.filename = fname
    sim.load()
    return(reconstruct(sim.root))

def save_h5(config, fname, update=False):

    cadetpath = run(['which', 'cadet-cli'], capture_output=True, text=True ).stdout.strip()
    Cadet.cadet_path = cadetpath

    sim = Cadet()
    sim.filename = fname

    if update:
        sim.root.update(Dict(config))
    else:
        sim.root = Dict(config)

    ## Do not add meta and output info to h5 file
    try:
        del(sim.root['meta'])
        del(sim.root['output'])
    except KeyError:
        pass

    sim.save()
    return sim

def delete_key(config:dict, key:str):
    try:
        del(config[key])
    except KeyError:
        pass
    return config


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("params", nargs='*', help="json/yaml parameter files")
    ap.add_argument("-m1", "--mono1d", action='store_true')
    ap.add_argument("-m2", "--mono2d", action='store_true')
    ap.add_argument("-p1", "--poly1d", action='store_true')
    ap.add_argument("-p2", "--poly2d", action='store_true')
    ap.add_argument("--mode", choices=['1d', '1d_dpfr', '2d_serial_inlet', '2d_parallel_inlet', '2d_serial_dpfr' ], required=True)
    args = vars(ap.parse_args())


    for paramfile in args['params']:
        
        paramfilepath = Path(paramfile)
        params = Dict(YAML(typ='safe').load(paramfilepath))

        outputbasefilename = paramfilepath.stem

        mode = params.mode
        if args['mode']:
            mode = args['mode']

        print("Loaded parameter file:", paramfile)

        if args['mono1d']:
            simulate(mode, 1, False, params, outputbasefilename + '.mono1d.h5')
        if args['poly1d']:
            simulate(mode, 1, True, params,  outputbasefilename + '.poly1d.h5')
        if args['mono2d']:
            simulate(mode, 2, False, params,  outputbasefilename + '.mono2d.h5')
        if args['poly2d']:
            simulate(mode, 2, True, params,  outputbasefilename + '.poly2d.h5')

if __name__ == "__main__":
    main()
