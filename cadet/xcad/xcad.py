#!/usr/bin/env python3

# TODO: WRITE OUT TO YAML
# TODO: mono2d
# TODO: Merge 2D and 1D GRM
# TODO: allow reading pickled files for volfracs, radii and porosities
# TODO: Make it smoother to use pack.py and this script
# TODO: Use nbeads and par_radius to calculate porosity in mono cases?
# TODO: Assumption is that one json file should represent the information for a full 3D column, If so, col_length/porosity with and without --no-dpfr will be different

"""
Run CADET model based on XNS model of chromatography.

unit_000: inlet
unit_001: (d)pfr
unit_002: GRM/GRM2D
unit_003: (d)pfr

Input file is in json format.

@example: `./xcad.py long.json -m1 -m2 -p1 -p2`

The above line runs the model with the parameters in long.json for all 4 variants: monodisperse 1d, monodisperse 2d, polydisperse 1d, polydisperse 2d

@NOTE: Information required for 2d/psd models can be generated via pack.py in the scripts dir. of this repo.
@NOTE: When running pack.py, consider that it might not be exactly similar to genmesh's output. For example, genmesh has the recently implemented porosity control to adjust column porosity for the polydisperse meshes. This isn't available in pack.py.

"""



from addict import Dict
import numpy as np
import argparse
from cadet import Cadet
from subprocess import run

import json
from ruamel.yaml import YAML
from pathlib import Path

# cadetpath = "/home/jayghoshter/local/bin/cadet-cli"

cadetpath = run(['which', 'cadet-cli'], capture_output=True, text=True ).stdout.strip()
Cadet.cadet_path = cadetpath


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

    ## Linear model
    GRMUnit.adsorption.lin_ka     = params.lin_ka
    GRMUnit.adsorption.lin_kd     = params.lin_kd

    ## MCL
    GRMUnit.adsorption.mcl_ka     = params.mcl_ka
    GRMUnit.adsorption.mcl_kd     = params.mcl_kd
    GRMUnit.adsorption.mcl_qmax   = params.mcl_qmax

    ## SMA
    GRMUnit.adsorption.sma_ka         = params.sma_ka
    GRMUnit.adsorption.sma_kd         = params.sma_kd
    GRMUnit.adsorption.sma_nu         = params.sma_nu
    GRMUnit.adsorption.sma_sigma      = params.sma_sigma
    GRMUnit.adsorption.sma_lambda     = params.sma_lambda
    GRMUnit.adsorption.sma_refc0      = params.sma_refc0
    GRMUnit.adsorption.sma_refq       = params.sma_refq

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

def create_simulation(mode, ndim, psd, params, outfile):

    model = Cadet()
    model.filename = outfile

    root                     = model.root

    column_unit=None
    if params.column_type == "GRM":
        if ndim == 2:
            column_unit = get_2dgrm_unit(params, psd=psd)
        else:
            column_unit = get_grm_unit(params, psd=psd)
    elif params.column_type == "LRM":
        column_unit = get_lrm_unit(params)

    # if ndim == 2:
    #     column_unit = get_2dgrm_unit(params, psd=psd)
    # else:
    #     column_unit = get_grm_unit(params, psd=psd)

    root.input['return'] = get_ioreturn()

    # Mode = "no_dpfr", "serial", "parallel 2d"
    if mode == "column":
        root.input.model.nunits  = 2
        root.input.model.unit_000 = get_inlet_unit(params)
        root.input.model.unit_001 = column_unit

        root.input['return'].unit_000 = get_ioreturn().unit
        root.input['return'].unit_001 = get_ioreturn().unit

        nports = [ root.input.model.get(unit).ports for unit in sorted(root.input.model.keys()) if 'unit_' in unit ]
        root.input.model.connections.switch_000.connections = getConnectionsMatrix(nports, params.velocity, params.radial_disc_type, params.col_radius)

    elif mode == "serial_dpfr":
        root.input.model.nunits   = 4
        root.input.model.unit_000 = get_inlet_unit(params)
        root.input.model.unit_001 = get_dpfr_unit(params)
        root.input.model.unit_002 = column_unit
        root.input.model.unit_003 = get_dpfr_unit(params)

        root.input['return'].unit_000 = get_ioreturn().unit
        root.input['return'].unit_001 = get_ioreturn().unit
        root.input['return'].unit_002 = get_ioreturn().unit
        root.input['return'].unit_003 = get_ioreturn().unit

        nports = [ root.input.model.get(unit).ports for unit in sorted(root.input.model.keys()) if 'unit_' in unit ]
        root.input.model.connections.switch_000.connections = getConnectionsMatrix(nports, params.velocity, params.radial_disc_type, params.col_radius)

    elif mode == "parallel_inlet":
        nrad     = len(params.col_porosity_2d)
        root.input.model.nunits   = nrad + 2
        for index in range(nrad):
            root.input.model[f"unit_{index:03d}"] = get_inlet_unit(params)
            root.input['return'][f"unit_{index:03d}"] = get_ioreturn().unit

        root.input.model[f"unit_{nrad:03d}"] = column_unit
        root.input['return'][f"unit_{nrad:03d}"] = get_ioreturn().unit

        root.input.model[f"unit_{nrad+1:03d}"] = get_outlet_unit(params)
        root.input['return'][f"unit_{nrad+1:03d}"] = get_ioreturn().unit

        nports = [ root.input.model.get(unit).ports for unit in sorted(root.input.model.keys()) if 'unit_' in unit ]
        # root.input.model.connections.switch_000.connections = getConnectionsParallel(max(nports), params.velocity, params.radial_disc_type, params.col_radius)

        connections = []
        conn1 = connectParallel([(x,0) for x in range(nrad)], [(nrad, y) for y in range(nrad)], params.velocity, params.radial_disc_type, params.col_radius)
        conn2 = connectSerial([(nrad, y) for y in range(nrad)], [(nrad+1,0)], params.velocity, params.radial_disc_type, params.col_radius) 
        connections.extend(conn1)
        connections.extend(conn2)
        root.input.model.connections.switch_000.connections = connections


    # nports = [ root.input.model.get(unit).ports for unit in sorted(root.input.model.keys()) if 'unit_' in unit ]

    root.input.model.connections.nswitches                 = 1
    root.input.model.connections.switch_000.section        = 0
    root.input.model.connections.connections_include_ports = 1

    # if "parallel" in mode:
    #     root.input.model.connections.switch_000.connections = getConnectionsParallel(max(nports), params.velocity, params.radial_disc_type, params.col_radius)
    # else:
    #     root.input.model.connections.switch_000.connections = getConnectionsMatrix(nports, params.velocity, params.radial_disc_type, params.col_radius)

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

# Assumes only serial connection, and indices starting from 0 onwards
def getConnectionsMatrix(nPorts, velocity, radial_disc_type, col_radius):
    ## To calc flowrates:
    ## > 1. rShells -> Area -> Flowrate for constant velocity input
    connections = []
    nUnits = len(nPorts)
    for i in range(nUnits-1):
        flowrates = getFlowrates(max(nPorts[i], nPorts[i+1]), velocity, radial_disc_type, col_radius)
        line = [ [i, i+1, x, y, -1, -1, flowrates[max(x,y)]] for x in range(nPorts[i]) for y in range(nPorts[i+1]) ]
        line = [x for sub in line for x in sub]
        connections.extend(line)

    return connections

def getFlowrates(nrad, velocity, radial_disc_type, col_radius):
    nRegions = nrad
    nShells = nRegions + 1 #Including r = 0
    rShells = []

    ## NOTE: move to current unit?
    # R = params.col_radius
    R = col_radius

    if radial_disc_type == 'EQUIVOLUME':
        for n in range(nShells):
            rShells.append(R * np.sqrt(n/nRegions))
    elif radial_disc_type == 'EQUIDISTANT':
        for n in range(nShells):
            rShells.append(R * (n/nRegions))

    flowrates = [velocity * np.pi *  (rShells[i+1]**2 - rShells[i]**2) for i in range(nRegions)]
    return flowrates

## Assumes indices go from (0 - n-1)  -> n, where n is the unit index of the multiport operation
def getConnectionsParallel(nrad, velocity, radial_disc_type, col_radius):
    connections = []
    # nSerialUnits = len(nPorts[0]) * len(nPorts[1])
    flowrates = getFlowrates(nrad, velocity, radial_disc_type, col_radius)
    for index in range(nrad):
        line = [ index, nrad, 0, index, -1, -1, flowrates[index] ]
        connections.extend(line)
    return connections

def connectParallel(unit_port_left:list, unit_port_right:list, velocity, radial_disc_type, col_radius):
    """ Connect two units in parallel
    unit_port_left and unit_port_right are lists of tuples: [(0,0), (1,0), (2,0)] and [(3,0), (3,1), (3,2)]
    """
    connections = []
    assert(len(unit_port_left) == len(unit_port_right))

    flowrates = getFlowrates(max(len(unit_port_left), len(unit_port_right)), velocity, radial_disc_type, col_radius) 

    for left, right, flowrate in zip(unit_port_left, unit_port_right, flowrates):
        line = [ left[0], right[0], left[1], right[1], -1, -1, flowrate]
        connections.extend(line)

    return connections

def connectSerial(unit_port_left:list, unit_port_right:list, velocity, radial_disc_type, col_radius):
    """ Connect two units serially 
    unit_port_left and unit_port_right are lists of tuples: [(3,0), (3,1), (3,2)] and [(4,0)]
    """
    connections = []

    # assert(sum([x[1] for x in unit_port_left]) == sum([x[1] for x in unit_port_right]))

    flowrates = getFlowrates(max(len(unit_port_left), len(unit_port_right)), velocity, radial_disc_type, col_radius) 

    pairs = [(x,y) for x in unit_port_left for y in unit_port_right]

    print(unit_port_left)
    print(unit_port_right)
    print(pairs)
    print(len(pairs))
    print(len(flowrates))

    assert(len(pairs) == len(flowrates))

    for pair, flowrate in zip(pairs, flowrates):
        line = [pair[0][0], pair[1][0], pair[0][1], pair[1][1], -1, -1, flowrate]
        connections.extend(line)

    return connections
        




    
def simulate(mode, ndim, psd, params, outfile):
    sim = create_simulation(mode, ndim, psd, params, outfile)
    sim.save()
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

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("params", nargs='*', help="json parameter files")
    ap.add_argument("-m1", "--mono1d", action='store_true')
    ap.add_argument("-m2", "--mono2d", action='store_true')
    ap.add_argument("-p1", "--poly1d", action='store_true')
    ap.add_argument("-p2", "--poly2d", action='store_true')
    ap.add_argument("--mode", choices=['column', 'parallel_inlet', 'serial_dpfr', 'parallel_dpfr'], required=True)
    args = vars(ap.parse_args())

    for paramfile in args['params']:

        params = Dict(YAML(typ='safe').load(Path(paramfile)))

        outputbasefilename = paramfile.replace('.json', '')

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
