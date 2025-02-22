from chromoo.cadetSimulation import CadetSimulation
import numpy as np
from joblib import Memory

from matplotlib import pyplot as plt
import matplotlib as mpl
from cycler import cycler

CACHE_DIR = f"cache"
memory = Memory(location=CACHE_DIR, verbose=0)

## Caching helps us iterate on post-proc faster. Delete CACHE_DIR directory for reset.
@memory.cache
def sim_load_run_post(unit_id):
    sim = CadetSimulation()
    sim.load_file('./simulation.yaml')
    sim.save()
    sim.run()
    sim.load()

    sim.post_mass_bulk(unit_id)
    # sim.post_mass_par(unit_id)
    # sim.post_mass_solid(unit_id)

    return sim

# mbulk_data = sim.root.output.post[unit_str].post_mass_bulk
# mpar_data = sim.root.output.post[unit_str].post_mass_par
# msolid_data = sim.root.output.post[unit_str].post_mass_solid

t,mbulk_ref = np.loadtxt('./reference/radial_shell_integrate_time_scalar_0_0_EQUIDISTANT_N1_CLIPPED_U.DV.csv', delimiter=',').T
t,mpar_ref = np.loadtxt('./reference/radial_shell_integrate_time_scalar_0_0_EQUIDISTANT_N1_FULL_U.DV_SCALED_0.25.csv', delimiter=',').T
t,msolid_ref = np.loadtxt('./reference/radial_shell_integrate_time_scalar_0_0_EQUIDISTANT_N1_FULL_U.DV_SCALED_0.75.csv', delimiter=',').T

def parameterized_runs(nsims):
    unit_id = 2
    unit_str = f'unit_{unit_id:03d}'

    ncols_1000 = list(range(10,1001, 100))
    ncols_void = list(range(1,101,10))

    for i, isim, in enumerate(range(nsims)):

        print(f"Running simulation #{isim:03d}")
        sim = CadetSimulation()
        sim.load_file('./simulation.yaml')

        # sim.root.input.model.unit_001.discretization.ncol = ncols_void[i]
        # sim.root.input.model.unit_003.discretization.ncol = ncols_void[i]
        # sim.root.input.model.unit_002.discretization.ncol = ncols_1000[i]
        sim.root.input.model.unit_002.col_dispersion = 4.765892089827942e-07
        sim.root.input.model.unit_002.film_diffusion = 8.099635961125705e-06

        sim.save()
        sim.run()
        sim.load()

        times = sim.root.input.solver.user_solution_times

        sim.post_mass_bulk(unit_id)
        mbulk_data = sim.root.output.post[unit_str].post_mass_bulk.sum(axis=1)

        sim.post_mass_par(unit_id)
        mpar_data = sim.root.output.post[unit_str].post_mass_par.sum(axis=1)

        sim.post_mass_solid(unit_id)
        msolid_data = sim.root.output.post[unit_str].post_mass_solid.sum(axis=1)

        ax.plot(times, mbulk_data + mpar_data + msolid_data, ls='dashed')

    return 


fig, ax = plt.subplots()
cmap = mpl.colormaps.get_cmap('tab20')
color_cycler = cycler(color=cmap.colors)
ax.set_prop_cycle(color_cycler)

ax.plot(t, mbulk_ref + mpar_ref + msolid_ref, color='black')

parameterized_runs(1) 
plt.legend()
plt.savefig('plot_run.pdf')

