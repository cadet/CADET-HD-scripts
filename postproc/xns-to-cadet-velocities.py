#!/usr/bin/env python3

"""xns-to-cadet-velocities.py

Due to 2DGRM using the compartment model, flow only happens along axial channels. 

This means that interstitial velocity in a given channel only depends on flowrate in that channel, not the entire inlet(s) to the full column. Thus v_{int} is inversely proportional to epsilon (column porosity)

Whereas it is known from the 3D simulations that the interstitial velocity is directly proportional to the column porosity. So in order to (correctly) run a cadet simulation in 2DGRM, we need to extrapolate the velocity profile from XNS into a corresponding INLET profile in CADET such that the interstitial velocities are the same.

THIS script is not required. XNS_radial_flowrates are equal to CADET_radial_flowrates.?
except that that's not the goal. 
XNS flowrates --> XNS velocities given XNS Areas/Porosities.
Double check that the final CADET velocities are the same as XNS velocities. Because the goal is to ensure the velocities are matching
"""


