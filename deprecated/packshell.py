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

# TODO: Read and center the column

import sys
import struct
import itertools
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
from math import asin
from scipy.special import ellipk, ellipe, ellipkinc, ellipeinc

def main():
    packing = sys.argv[1]
    zBot = float(sys.argv[2])
    zTop = float(sys.argv[3])
    scaling_factor = float(sys.argv[4])
    # rFactor = 0.9997
    rFactor = 1
    meshScalingFactor = 1e-4

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

    r = [item/2 * scaling_factor * rFactor * meshScalingFactor for item in d]
    x = [item * scaling_factor * meshScalingFactor for item in x]
    y = [item * scaling_factor * meshScalingFactor for item in y]
    z = [item * scaling_factor * meshScalingFactor for item in z]

    x, y = centerCylinder(x,y,z,r)
    xmax, xmin, ymax, ymin, zmax, zmin = updateBounds(x,y,z,r)
    R = max(xmax, -xmin, ymax, -ymin)
    rlimit = R/np.sqrt(2)

    # NOTE: Bed length is not similar to that calculated in genmesh
    h = zmax - zmin

    inner = []
    outer = []
    onbound = []
    inbound = []
    outbound = []
    allbound = []
    for i in range(len(r)):
        # if ((x[i]**2 + y[i]**2)  <= (R**2)/2):
        if (np.sqrt(x[i]**2 + y[i]**2) + r[i]  < rlimit):
            inner.append(r[i])
        elif (np.sqrt(x[i]**2 + y[i]**2) - r[i]  > rlimit):
            outer.append(r[i])
        else:
            if ((x[i]**2 + y[i]**2)  < rlimit**2):
                inbound.append(r[i])
            elif ((x[i]**2 + y[i]**2)  > rlimit**2):
                outbound.append(r[i])
            elif ((x[i]**2 + y[i]**2)  == rlimit**2):
                onbound.append(r[i])
            allbound.append(r[i])

    print(len(inner), len(outer), len(inbound), len(outbound), len(onbound))
    print(sum([len(inner), len(outer), len(inbound), len(outbound), len(onbound)]))

    # print("INNER")
    # histo(inner, 'inner')
    # print("OUTER")
    # histo(outer, 'outer')

    Vol_inner = sum([4/3 * np.pi * r**3 for r in inner])
    Vol_outer = sum([4/3 * np.pi * r**3 for r in outer])

    # print(Vol_inner, Vol_outer)
    # sys.exit()


    for i in range(len(allbound)):
        rad = allbound[i]
    # for i in range(len(outbound)):
    #     rad = outbound[i]
    # for i in range(len(inbound)):
    #     rad = inbound[i]
        tempRho = rlimit/rad
        print(tempRho)
        tempEta = np.sqrt(x[i]**2 + y[i]**2)/rad
        print(tempRho, tempEta, tempEta - tempRho)
        Vint = CylSphIntVolume(tempRho, tempEta)
        Vint = Vint * rad**3
        Vsph = 4/3 * np.pi * rad**3
        Vol_inner = Vol_inner + Vint
        Vol_outer = Vol_outer + (Vsph - Vint)

    print(Vol_inner, Vol_outer)
    VolInCyl = np.pi * rlimit**2 * h
    VolOutCyl = np.pi * (R**2 - rlimit**2) * h
    print(Vol_inner/VolInCyl, Vol_outer/VolOutCyl)


    # TODO: scale cylsphvolume for sphere radius and cylinder radius
    # TODO: MS OUT and MS IN cases separately
    # TODO: Implement

def centerCylinder(x, y, z, r):
    xmax, xmin, ymax, ymin, _, _ = updateBounds(x, y, z, r)
    offsetx = -(xmax + xmin)/2
    offsety = -(ymax + ymin)/2
    x = [ xc + offsetx for xc in x ]
    y = [ yc + offsety for yc in y ]
    return x, y

def updateBounds(x, y, z, r):
    xpr = []
    xmr = []
    ypr = []
    ymr = []
    zpr = []
    zmr = []
    for i in range(len(x)):
        xpr.append(x[i]+r[i])
        xmr.append(x[i]-r[i])
        ypr.append(y[i]+r[i])
        ymr.append(y[i]-r[i])
        zpr.append(z[i]+r[i])
        zmr.append(z[i]-r[i])

    xmax = max(xpr)
    ymax = max(ypr)
    ymin = min(ymr)
    xmin = min(xmr)
    zmax = max(zpr)
    zmin = min(zpr)

    return xmax, xmin, ymax, ymin, zmax, zmin

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


def histo(radii, filename):
    V=[4*np.pi*x*x*x/3 for x in radii]
    h,e = np.histogram(radii, bins=20, density=True, weights=V)

    frac=[x/sum(h) for x in h]
    print(sum(frac))
    print(frac)
    # print(h)
    # print(e)
    # print(np.diff(e))
    w=2
    avg=np.convolve(e, np.ones(w), 'valid') / w
    # print(list(avg))


    with plt.style.context(['science']):
        matplotlib.rcParams['font.sans-serif'] = "Verdana"
        matplotlib.rcParams['font.family'] = "sans-serif"

        # # plt.rc('axes', prop_cycle=(cycler('color', ['#0C5DA5', '#00B945', '#FF9500', '#FF2C00', '#845B97', '#474747', '#9e9e9e']) ))
        # plt.rc('axes', prop_cycle=(cycler('color', ['#0C5DA5', '#00B945']) ))

        fig, ax = plt.subplots()
        ax.hist(radii, bins=20)

        # ax.legend(loc='lower center', bbox_to_anchor=(0.5,-0.7), shadow=True)
        ax.set(title=filename)
        # ax.set(xlabel='Time (s)', fontname="Arial")
        ax.set_xlabel('Bead Radius ($m$)')
        ax.set(ylabel='Frequency')
        # ax.autoscale(tight=True)
        fig.savefig(filename + '.pdf')
        # fig.savefig(filename+'.jpg', dpi=300)

def CylSphIntVolume(rho, eta):
    nu = asin(eta - rho)
    k = np.sqrt((1- (eta - rho)**2)/(4*rho*eta))
    V = 2/3 * np.pi * ( 1 - L0(nu, k)) - 8/9 * np.sqrt(rho * eta) * (6 * rho**2 + 2 * rho * eta - 3) * (1 - k**2) * ellipk(k)  + 8/9 * np.sqrt(rho * eta) * (7 * rho**2 + eta**2 - 4) * ellipe(k)
    return V

def L0(xi, k):
    kdash = np.sqrt(1 - k**2)
    # print(ellipe(k), ellipkinc(xi, kdash), ellipk(k),)
    value = 2/np.pi * (ellipe(k) * ellipkinc(xi, kdash) + ellipk(k) * ellipeinc(xi, kdash) - ellipk(k) * ellipkinc(xi, kdash))
    return value

if __name__ == "__main__":
    print(sys.version)
    print(__doc__)
    main()
