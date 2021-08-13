#!/usr/bin/python3

"""
Calculate lateral porosity profile for periodic packings.

Prerequisities:
    - overlap: `pip install overlap` (See https://github.com/severinstrobl/overlap)
    - numpy

Usage:
    ./script.py <beads.xyzd> -n <num x-slices> -b <box dimensions in x0 y0 z0 dx dy dz>

    |width of slice, advancing -->
    |--|
--------------------
|   .  .           |
|   .  .           |
|   .  .           |
|   .  .           |
--------------------

|-------dX----------|

"""

import struct
import numpy as np
import argparse
import overlap
from dataclasses import dataclass
from itertools import tee, islice

@dataclass(init=True, order=True, repr=True, frozen=True)
class Sphere:
    x: float
    y: float
    z: float
    r: float

def bin_to_arr(filename, format):
    """
    Read binary data into array
    """

    with(open(filename, 'rb')) as input:
        myiter = struct.iter_unpack(format, input.read())

        arr = []
        for i in myiter:
            arr.append(i[0])

        return arr

def grouper(iterable, n):
    """
    Group binary data into chunks after reading
    """
    it = iter(iterable)
    while True:
       chunk = tuple(islice(it, n))
       if not chunk:
           return
       yield chunk

class PackedBed:
    def __init__(self, filename, ):
        self.fname = filename
        self.dataformat = '<d'
        self.read_packing()

    def read_packing(self):
        """
        Read a packing from a little endian double-sized binary file

        The input contains "duplicated" beads, those that peek out of one periodic boundary and into the other
        The overlap algorithm takes care of not accounting for particle volume peeking out of the periodic boundaries
        """
        self.beads = []
        arr = bin_to_arr(self.fname, self.dataformat)
        for index, chunk in enumerate(grouper(arr,4)):
            x = chunk[0]
            y = chunk[1]
            z = chunk[2]
            r = chunk[3]/2
            self.beads.append(Sphere(x, y, z, r))

@dataclass(init=True, order=True, repr=True, frozen=True)
class Hexahedron:
    vertices: tuple


def hex_sph_x(hex, sph):
    _hex = overlap.Hexahedron(hex.vertices)

    ## ((x,y,z), r)
    _sphere = overlap.Sphere((sph.x, sph.y, sph.z), sph.r)

    # full_sphere_vol = 4/3 * np.pi * (sph.r**3)
    # print(f"{full_sphere_vol = }")
    # print(f"{full_sphere_vol/2 = }")
    # print(f"{full_sphere_vol/4 = }")
    # print(f"{full_sphere_vol/8 = }")
    result = overlap.overlap(_sphere, _hex)
    return result


def pairwise(iterable):
    """Iterate in pairs

    >>> list(pairwise([0, 1, 2, 3]))
    [(0, 1), (1, 2), (2, 3)]
    >>> tuple(pairwise([])) == tuple(pairwise('x')) == ()
    True
    """
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)

def csvWriter(filename, x, y):
    import csv
    with open(filename, 'w') as f:
        writer = csv.writer(f)
        writer.writerows(zip(x, y))

def main():

    ap = argparse.ArgumentParser()
    ap.add_argument('-b', '--box', nargs=6, type=float, required=True, help='x0 y0 z0 dx dy dz for the box')
    ap.add_argument('-n', '--nRegions', default=10, type=int, help='Number of hex slices of the box to calculate over')
    ap.add_argument('beadsfile', help="Beads file in xyzd little endian binary")
    args = vars(ap.parse_args())

    X0, Y0, Z0, dX, dY, dZ = args['box']

    X1 = X0+dX
    Y1 = Y0+dY
    Z1 = Z0+dZ

    nregions = args['nRegions']
    ngrid = nregions + 1
    xgrid = np.linspace(X0, X1, ngrid)

    packedBed = PackedBed(args['beadsfile'])

## x-direction, 0 -> end
    xs = []
    ys = []
    for xl, xr in pairwise(xgrid):
        hex = Hexahedron((
            (xl,Y0,Z0),
            (xr,Y0,Z0),
            (xr,Y1,Z0),
            (xl,Y1,Z0),
            (xl,Y0,Z1),
            (xr,Y0,Z1),
            (xr,Y1,Z1),
            (xl,Y1,Z1),
            ))

        slice_sph_vol = sum((hex_sph_x(hex, sphere) for sphere in packedBed.beads))
        slice_hex_vol = dY * dZ * (xr - xl)
        slice_porosity = (1 - slice_sph_vol/slice_hex_vol)

        xs.append((xl+xr)/2)
        ys.append(slice_porosity)

    csvWriter('porosities.csv', xs, ys)

if __name__ == "__main__":
    main()
