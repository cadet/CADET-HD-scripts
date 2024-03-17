#!/usr/bin/env python3

"""
Script to convert outputs of paravision radial_shell_integrate into cadet form. 
In radial_shell_integrate, we typically output t,y data per radial zone (as separate files) with the shape of (t). 
Cadet stores solutions with time as first axis, (t, n_ax, n_rad, n_sol...) etc. 

Doing this allows defining all 5 objectives as a super-objective in the chromoo config:

```chromoo.yaml
objectives: 
    - name: $M_{b}$
      filename: ./reference/bulk/y.out
      times: ./reference/bulk/x.out
      score: nrmse
      path: output.post.unit_001.post_mass_bulk
      sum_data_axis: 1
      shape: [99,5]
```
"""

import numpy as np
from pathlib import Path

files = list(sorted(Path('.').glob('./radial_shell_integrate_time_scalar_0_*_GRM2D_FULL_U.DV.csv')))

ms = []
ts = []
for f in files:
    t,m = np.loadtxt(f,delimiter=',').T
    ms.append(m)
    ts.append(t)

ms_as_cadet = np.array(ms).T
np.savetxt('y.out', ms_as_cadet.ravel())
np.savetxt('x.out', ts[0])
