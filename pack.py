#!/usr/bin/env python3

"""
@name: pack.py
@desc: Read binary packing.xyzd files, and generate porosity profiles.
@theory: Cylinder/Sphere intersection volume equations from http://dx.doi.org/10.1016/s1385-7258(61)50049-2
@usage: ./pack.py <packing.xyzd> <zBot after scaling> <zTop after scale> <scaling factor>
        zBot: bottom limit of slice to look for beads by center point.
        zTop: top limit of slice ...
        scaling factor: two packing.xyzd might not have the same dimensions. This helps fix that.

"""

# DONE: implement use of histo
# INPROGRESS: Documentation
# DONE: Generate results for CADET input: binned radii & volume fractions
# TODO: argparse, argument handling
# TODO: read genmesh input file??
# TODO: Handle rho == eta edge cases
# TODO: Auto handle scaling_factor: (Updatebounds, scale to fit Cyl Radius = 5)
# TODO: Better parallelization?

import sys
import struct
import itertools
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
from math import asin,sqrt,pi
from mpmath import ellipk, ellipe, ellipf, nstr
from multiprocessing import Pool
from functools import partial

class Bead:
    """Class for individual beads"""

    def __init__(self, x, y, z, r):
        self.x = x
        self.y = y
        self.z = z
        self.r = r

    def pos(self):
        return np.sqrt(self.x**2 + self.y**2)

    def volume(self):
        return 4/3 * np.pi * self.r**3

    def distance(self, other):
        return sqrt((self.x-other.x)**2 + (self.y-other.y)**2 + (self.z-other.z)**2)



class PackedBed:
    """Class for packed bed of beads. Can apply transformations on beads"""

    def __init__(self):
        self.beads = []

    def add(self, bead):
        self.beads.append(bead)

    def size(self):
        return len(self.beads)

    def volume(self):
        vol = 0
        for bead in self.beads:
            vol = vol + bead.volume()
        return vol

    def updateBounds(self):
        """
        Calculate bounding points for the packed bed.
        """

        xpr = []
        xmr = []
        ypr = []
        ymr = []
        zpr = []
        zmr = []

        for bead in self.beads:
            xpr.append(bead.x + bead.r)
            xmr.append(bead.x - bead.r)
            ypr.append(bead.y + bead.r)
            ymr.append(bead.y - bead.r)
            zpr.append(bead.z + bead.r)
            zmr.append(bead.z - bead.r)

        self.xmax = max(xpr)
        self.ymax = max(ypr)
        self.ymin = min(ymr)
        self.xmin = min(xmr)
        self.zmax = max(zpr)
        self.zmin = min(zpr)
        self.R = max(self.xmax, -self.xmin, self.ymax, -self.ymin)
        self.h = self.zmax - self.zmin
        self.CylinderVolume = pi * self.R**2 * self.h

    def moveBedtoCenter(self):
        """
        Translate bed center to origin of coordinate system.
        """
        self.updateBounds()
        offsetx = -(self.xmax + self.xmin)/2
        offsety = -(self.ymax + self.ymin)/2
        for bead in self.beads:
            bead.x = bead.x + offsetx
            bead.y = bead.y + offsety
        self.updateBounds()


def grouper(iterable, n):
    """Group binary data into chunks after reading"""
    it = iter(iterable)
    while True:
       chunk = tuple(itertools.islice(it, n))
       if not chunk:
           return
       yield chunk

def bin_to_arr(filename, f):
    """Read binary data into array"""

    with(open(filename, 'rb')) as input:
        myiter = struct.iter_unpack(f, input.read())

        arr = []
        for i in myiter:
            arr.append(i[0])

        return arr

def histo(radii, filename):
    """Create histogram for a particular bead size distribution.
        Also output volume fractions & mean radii to be used in CADET Polydisperse"""
    V=[4*np.pi*x*x*x/3 for x in radii]
    h,e = np.histogram(radii, bins=20, density=True, weights=V)

    frac=[x/sum(h) for x in h]
    # print(sum(frac))
    w=2
    avg=np.convolve(e, np.ones(w), 'valid') / w
    print('vol_frac:', frac)
    print('mean_radii:', list(avg))


    with plt.style.context(['science']):
        matplotlib.rcParams['font.sans-serif'] = "Verdana"
        matplotlib.rcParams['font.family'] = "sans-serif"

        fig, ax = plt.subplots()
        ax.hist(radii, bins=20)

        ax.set(title=filename)
        ax.set_xlabel('Bead Radius ($m$)')
        ax.set(ylabel='Frequency')
        fig.savefig(filename + '.pdf')

def CylSphIntVolume(rho, eta):
    """ Analytical Formulae to calculate intersection between cylinder and sphere.
        See http://dx.doi.org/10.1016/s1385-7258(61)50049-2 for more info.
    """
    if rho == 0.0:
        return 0
    elif (eta - rho) <= -1:
        return 4/3 * pi
    elif (eta - rho) >= 1:
        return 0

    ## NOTE: Ideally eta & rho are floats & never equal. But test cases are not handled yet. Similarly rho+eta == 1
    if eta == rho:
        print("Rho & Eta are Equal")

    if eta == 0 and 0 <= rho <= 1:
        V = 4/3 * pi - 4/3 * pi * (1 - rho**2)**(3/2)
        return V
    elif (rho + eta > 1):
        nu = asin(eta - rho)
        m = (1-(eta - rho)**2)/(4*rho*eta)

        K = ellipk(m)
        E = ellipe(m)

        F = ellipf(nu ,1-m)
        Ep = ellipe(nu, 1-m)

        L0 = 2/pi * (E * F + K * Ep - K * F )

        # V = (2/3 * pi * ( 1 - L0(nu, m) ) )\
        V = (2/3 * pi * ( 1 - L0 ) )\
        - (8/9 * sqrt(rho * eta) * (6 * rho**2 + 2 * rho * eta - 3) * (1 - m) * K)\
        + (8/9 * sqrt(rho * eta) * (7 * rho**2 + eta**2 - 4) * E)

        return V

    elif (rho + eta < 1):
        nu = asin((eta - rho)/(eta + rho))
        m = 4*rho*eta / (1 - (eta-rho)**2)
        K = ellipk(m)
        E = ellipe(m)
        F = ellipf(nu ,1-m)
        Ep = ellipe(nu, 1-m)
        L0 = 2/pi * (E * F + K * Ep - K * F )

        V = (2/3 * pi * ( 1 - L0 ))\
        - (4 * sqrt(1 - (eta-rho)**2) / (9*(eta+rho)) ) * (2*rho - 4*eta + (eta+rho)*(eta-rho)**2) * (1-m) * K\
        + (4/9 * sqrt(1 - (eta-rho)**2) * (7*rho**2 + eta**2 - 4) * E)

        return V

    else:
        print("ERROR")
        return 0

def main():
    packing = sys.argv[1]
    zBot = float(sys.argv[2])
    zTop = float(sys.argv[3])
    scaling_factor = float(sys.argv[4])
    rFactor = 0.99
    print("rFactor: ", rFactor)
    # rFactor = 1
    meshScalingFactor = 1e-4
    relativeBridgeRadius = 0.2

    fullBed = PackedBed()

    dataformat = "<f"
    arr = bin_to_arr(packing, dataformat)
    for chunk in grouper(arr,4):
        if (chunk[2] >= zBot/scaling_factor) and (chunk[2] <= zTop/scaling_factor):
            x = chunk[0] * scaling_factor * meshScalingFactor
            y = chunk[1] * scaling_factor * meshScalingFactor
            z = chunk[2] * scaling_factor * meshScalingFactor
            r = chunk[3]/2 * scaling_factor * rFactor * meshScalingFactor
            fullBed.add(Bead(x, y, z, r))

    fullBed.moveBedtoCenter()
    # R = fullBed.R + 0.01*meshScalingFactor ## Adding Rcyldelta
    R = fullBed.R
    h = fullBed.h

    print("R:", R)
    print("h:", h)
    print("Cylinder Volume:", fullBed.CylinderVolume)
    print("Packed Bed Volume:", fullBed.volume())
    print("nBeads: ", len(fullBed.beads))

    # bridgeOffsetRatio = 0.95
    bridgeOffsetRatio = sqrt(1 - relativeBridgeRadius**2)
    bridgeTol = 0.04 * meshScalingFactor

    addedBridgeVol, removedBridgeVol = bridgeVolumes(fullBed.beads, bridgeTol, relativeBridgeRadius, bridgeOffsetRatio)
    print("Bridge Added Volume:", addedBridgeVol)
    print("Capped Removed Volume:", removedBridgeVol)

    sys.exit(0)

    histo([bead.r for bead in fullBed.beads], 'Full Bed')

    nRegions = 100
    nShells = nRegions + 1 #Including r = 0
    rShells = []

    # shellType = 'equiVolume'
    shellType = 'equiRadius'

    if shellType == 'equiVolume':
        for n in range(nShells):
            rShells.append(R * sqrt(n/nRegions))
    elif shellType == 'equiRadius':
        for n in range(nShells):
            rShells.append(R * (n/nRegions))

    print("rShells:", rShells)

    volRegions = [0] * nRegions

    ## Multiprocessing code.
    ##      Create a partial function of volShellRegion(beads, rShells, i) --> parfunc(i)
    ##      map each 'i' to each process
    pool = Pool()
    parfunc = partial(volShellRegion, fullBed.beads, rShells)
    volRegions = pool.map(parfunc, range(nRegions))
    pool.close()
    pool.join()

    # print(volShellRegion(fullBed.beads, rShells, 0))
    volRegions = [float(item) for item in volRegions]

    volCylRegions = [pi * h * (rShells[i+1]**2 - rShells[i]**2) for i in range(nRegions)]
    porosities = [ float(1-n/m) for n,m in zip(volRegions, volCylRegions) ]
    avg_radius = [ (rShells[i] + rShells[i+1])/2 for i in range(nRegions) ]

    for i in range(nRegions):
        print(avg_radius[i], volRegions[i], volCylRegions[i], porosities[i])

    plotter(avg_radius, porosities, '', 'plot.pdf')


def volShellRegion(beads, rShells, i):
    """
    Find the intersection volumes between rShells[i] & rShells[i+1]
    """
    volShell=0
    for bead in beads:
        volBead = volBeadSlice(bead, rShells[i], rShells[i+1])
        volShell = volShell + volBead
    return volShell

def bridgeVolumes(beads, bridgeTol, relativeBridgeRadius, bridgeOffsetRatio):
    addedBridgeVol = 0
    removedBridgeVol = 0
    count = 0
    beadsCopy = beads.copy()
    for bead1 in beads:
        beadsCopy.remove(bead1)
        for bead2 in beadsCopy:
            beadDistance = bead1.distance(bead2)
            if beadDistance < bead1.r + bead2.r + bridgeTol:
                count = count + 1
                bridgeRadius = relativeBridgeRadius * min(bead1.r, bead2.r)
                intVol1 = volBridgeSlice(bead1, bridgeRadius, bridgeOffsetRatio)
                intVol2 = volBridgeSlice(bead2, bridgeRadius, bridgeOffsetRatio)
                addedBridgeVol = addedBridgeVol + pi * bridgeRadius**2 * (beadDistance - bridgeOffsetRatio * bead1.r - bridgeOffsetRatio * bead2.r) - intVol1 - intVol2
                removedBridgeVol = removedBridgeVol + intVol1 + intVol2
                ## NOTE: Some beads will be intersecting due to single precision. That's not handled here.
    print("Number of Bridges:", count)
    return addedBridgeVol, removedBridgeVol

def volBridgeSlice(bead, bridgeRadius, offsetRatio):
    """
    Volume of intersection between bridge and bead
    """
    rho = bridgeRadius/bead.r
    # eta = bead.pos()/bead.r ##FIXME, eta == 0
    eta = 0
    vol = CylSphIntVolume(rho, eta) * bead.r**3
    ## There's no need to find the accurate internal union volume since it will be deleted to find only the extra volume added by bridges in the first place.
    vol = vol/2 - pi * bridgeRadius**2 * offsetRatio * bead.r
    return vol

def volBeadSlice(bead, rInnerShell, rOuterShell):
    """
    Find intersection volume of an individual bead between two shells (cylinders)
    """
    rhoOuter = rOuterShell/bead.r
    etaOuter = bead.pos()/bead.r
    volOuter = CylSphIntVolume(rhoOuter, etaOuter) * bead.r**3
    rhoInner = rInnerShell/bead.r
    etaInner = bead.pos()/bead.r
    volInner = CylSphIntVolume(rhoInner, etaInner) * bead.r**3
    volIntBead = volOuter - volInner
    return volIntBead

def plotter(x, y, title, filename):
    with plt.style.context(['science']):
        fig, ax = plt.subplots()
        ax.plot(x, y)
        # legend = ax.legend(loc='best', shadow=True)
        ax.set(title=title)
        ax.set(xlabel='Radius')
        ax.set(ylabel='Porosity')
        ax.autoscale(tight=True)
        ax.set_ylim(0,1)
        fig.savefig(filename)


if __name__ == "__main__":
    print(sys.version)
    print(__doc__)
    main()
