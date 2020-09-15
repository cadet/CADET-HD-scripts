#!/usr/bin/env python3

"""
@name: packbin.py
@desc: create histogram bins from packing
@usage: ./packbin.py <packing.xyzd> zBot zTop preScalingFactor
@example: ./packbin.py poly-full.xyzd 0 15.10 2.1244954
@warning: === EXPERIMENTAL! ===
"""

#TODO: standardize
#TODO: clean output
#TODO: clean inputs, use argparse

import sys
import struct
import itertools
import numpy as np
from matplotlib import pyplot as plt
import matplotlib

def main():
    packing = sys.argv[1]
    zBot = float(sys.argv[2])
    zTop = float(sys.argv[3])
    scaling_factor = float(sys.argv[4])
    rFactor = 0.9997
    meshScalingFactor = 1e-4

    # print(zBot)
    # print(zTop)
    # print(scaling_factor)

    dataformat = "<f"
    arr = bin_to_arr(packing, dataformat)
    # x , y, z ,d = numpy.array()
    x = []
    y = []
    z = []
    d = []
    for chunk in grouper(arr,4):
        # print("\t".join("%.6E" % x for x in chunk))
        if (chunk[2] >= zBot/scaling_factor) and (chunk[2] <= zTop/scaling_factor):
            x.append(chunk[0])
            y.append(chunk[1])
            z.append(chunk[2])
            d.append(chunk[3])


    r = [x/2 * scaling_factor * rFactor * meshScalingFactor for x in d]
    V=[4*np.pi*x*x*x/3 for x in r]
    h,e = np.histogram(r, bins=10, density=True, weights=V)

    frac=[x/sum(h) for x in h]
    # print(sum(frac))
    print(frac)
    # print(h)
    # print(e)
    # print(np.diff(e))
    w=2
    avg=np.convolve(e, np.ones(w), 'valid') / w
    print(list(avg))


    # with plt.style.context(['science']):
    #     matplotlib.rcParams['font.sans-serif'] = "Verdana"
    #     matplotlib.rcParams['font.family'] = "sans-serif"

    #     # # plt.rc('axes', prop_cycle=(cycler('color', ['#0C5DA5', '#00B945', '#FF9500', '#FF2C00', '#845B97', '#474747', '#9e9e9e']) ))
    #     # plt.rc('axes', prop_cycle=(cycler('color', ['#0C5DA5', '#00B945']) ))

    #     fig, ax = plt.subplots()
    #     ax.hist(r, bins=10)

    #     # ax.legend(loc='lower center', bbox_to_anchor=(0.5,-0.7), shadow=True)
    #     # ax.set(title='Chromatogram')
    #     # ax.set(xlabel='Time (s)', fontname="Arial")
    #     ax.set_xlabel('Bead Radius ($m$)')
    #     ax.set(ylabel='Frequency')
    #     # ax.autoscale(tight=True)
    #     fig.savefig('histogram.pdf')
    #     # fig.savefig(filename+'.jpg', dpi=300)

    # V_beads=0
    # for i in range(len(d)):
    #     r = d[i]/2
    #     V_beads+=(4.0/3.0)*numpy.pi * r * r * r
    # xMax = max(x)
    # xMin = min(x)
    # yMax = max(y)
    # yMin = min(y)
    # zMax = max(z)
    # zMin = min(z)
    # print("zMax:", zMax)
    # print("zMin:", zMin)
    # print("xRad", (xMax-xMin)/2)
    # print("yRad", (yMax-yMin)/2)
    # R = (xMax-xMin)/2
    # l_bed = zMax-zMin
    # V_bed = numpy.pi * R * R * l_bed
    # porosity = (V_bed-V_beads)/V_bed
    # # print(V_beads)
    # # print(V_bed)
    # print("nbeads:",len(d))
    # print("porosity:", porosity)


def grouper(iterable, n):
    it = iter(iterable)
    while True:
       chunk = tuple(itertools.islice(it, n))
       if not chunk:
           return
       yield chunk

def bin_to_arr(filename, f):
    with(open(filename, 'rb')) as input:
        myiter = struct.iter_unpack(f, input.read())

        arr = []
        for i in myiter:
            arr.append(i[0])

        return arr



if __name__ == "__main__":
    import sys
    print(sys.version)
    print(__doc__)
    main()
