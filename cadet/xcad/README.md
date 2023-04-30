# Cadet XNS setup

**DEPRECATED: PLEASE USE GENCADET INSTEAD**

Cadet Python scripts for simulating ChromaHD models for comparison with XNS.

- [DONE] Improve interface
- [DONE] JSON input files?
- [TASK] Implement mono 2d

# Model
The Cadet model consists of 4 units: 
    - Inlet
    - DPFR
    - GRM/2DGRM
    - DPFR

# Usage

Use `pack.py` to generate all the porosity/vol_frac info required.
The script (xcad.py) uses a json input file to specify model parameters (lengths, porosities, binding parameters etc). 

This line should run all variants of the model: mono/poly x 1d/2d
`./xcad.py input.json -m1 -m2 -p1 -p2`

The `--no-dpfr` flag removes the dpfr units.

