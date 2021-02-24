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
# DONE: argparse, argument handling
# DONE: read genmesh input file??
# TODO: Handle rho == eta edge cases
# TODO: Auto handle scaling_factor: (Updatebounds, scale to fit Cyl Radius = 5)
# TODO: Better parallelization?

import sys
import struct
import itertools
import numpy as np
import argparse
from matplotlib import pyplot as plt
import matplotlib
from math import asin,sqrt,pi
from mpmath import ellipk, ellipe, ellipf, nstr
from multiprocessing import Pool
from functools import partial
import pickle
import os.path


# NBINS    = 20                   ## NPARTYPE in CADET
# NREGIONS = 100                   ## NRAD in CADET

NBINS    = 10                   ## NPARTYPE in CADET
NREGIONS = 20                   ## NRAD in CADET

def csvWriter(filename, x, y):
    import csv
    with open(filename, 'w') as f:
        writer = csv.writer(f)
        writer.writerows(zip(x, y))

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

        radList = [ bead.r for bead in self.beads ]
        self.rmax = max(radList)
        self.rmin = min(radList)
        self.ravg = sum(radList)/len(radList)

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

def pickler(item, filename):
    """
    > Save item to file.
    > For easy transport of lists/arrays to cadet
    """
    with open(filename, "wb") as f:
        pickle.dump(item, f)
        print("pickled", filename)

def unpickler(filename):
    with open(filename, "rb") as f:
        item = pickle.load(f)
        print("Unpickled", filename)
    return item

def histo(radii, **kwargs):
    """Create histogram for a particular bead size distribution.
        Also output volume fractions & mean radii to be used in CADET Polydisperse"""

    filename = kwargs.get('filename', None)
    bins = kwargs.get('bins', NBINS)

    V=[4*np.pi*x*x*x/3 for x in radii]
    # print(radii)
    h,e = np.histogram(radii, bins=bins, density=True, weights=V)
    # print("H: ", h)
    frac=[x/sum(h) for x in h]
    # print(sum(frac))
    w=2
    avg=np.convolve(e, np.ones(w), 'valid') / w

    if filename:
        with plt.style.context(['science']):
            matplotlib.rcParams['font.sans-serif'] = "Verdana"
            matplotlib.rcParams['font.family'] = "sans-serif"

            fig, ax = plt.subplots()
            ax.hist(radii, bins=bins)

            ax.set(title=filename)
            ax.set_xlabel('Bead Radius ($m$)')
            ax.set(ylabel='Frequency')
            fig.savefig(filename + '.pdf')

    return frac, list(avg)

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

# def getInput(infile, key):
#     pass

def main():

    ap = argparse.ArgumentParser()

    ap.add_argument("-i", "--inputfile", help="Input file used for genmesh")
    ap.add_argument("-p", "--packing", help="Packing file")
    ap.add_argument("-z", "--zlimits", type=float, nargs=2, help="bottom z value for packing slice")
    ap.add_argument("-psf", "--pre-scaling-factor", type=float, default=1, help="Scaling factor to account for differences in mono/poly packing data")
    ap.add_argument("-msf", "--mesh-scaling-factor", type=float, default=1e-4, help="Post meshing scaling factor")
    ap.add_argument("-rf", "--r-factor", type=float, default=1, help="Bead radius shrinking factor")

    # ap.add_argument("FILES", nargs='*', help="files..")

    args = vars(ap.parse_args())

    infiledict = {}

    if args['inputfile']:
        fp = open(args['inputfile'], 'r')
        infiledict = {line.strip().split()[0].strip(): line.strip().split()[1].strip() for line in fp}
        fp.close()

        packing = infiledict['packing']
        zBot = float(infiledict['zBot'])
        zTop = float(infiledict['zTop'])
        scaling_factor = float(infiledict['preScalingFactor'])
        rFactor = float(infiledict['rFactor'])
        meshScalingFactor = float(infiledict['Mesh.ScalingFactor'])
    else:
        packing = args['packing']
        zBot = args['zlimits'][0]
        zTop = args['zlimits'][1]
        scaling_factor = args['pre_scaling_factor']
        rFactor = args['r_factor']
        meshScalingFactor = args['mesh_scaling_factor']

    print("packing:", packing)
    print("zBot:", zBot)
    print("zTop:", zTop)
    print("scaling_factor:", scaling_factor)
    print("rFactor:", rFactor)
    print("meshScalingFactor:", meshScalingFactor)

    relativeBridgeRadius = 0.2

    rCylDelta = 0.01*meshScalingFactor
    inlet = 2.5 * meshScalingFactor
    outlet = 2.5 * meshScalingFactor

    fullBed = PackedBed()

    dataformat = "<f" ## For old packings with little endian floating point data. Use <d for new ones
    arr = bin_to_arr(packing, dataformat)
    for chunk in grouper(arr,4):
        if (chunk[2] >= zBot/scaling_factor) and (chunk[2] <= zTop/scaling_factor):
            x = chunk[0] * scaling_factor * meshScalingFactor
            y = chunk[1] * scaling_factor * meshScalingFactor
            z = chunk[2] * scaling_factor * meshScalingFactor
            r = chunk[3]/2 * scaling_factor * rFactor * meshScalingFactor
            fullBed.add(Bead(x, y, z, r))

    fullBed.moveBedtoCenter()
    R = fullBed.R + rCylDelta ## Adding Rcyldelta
    h = (zTop - zBot)*meshScalingFactor  + inlet + outlet
    hBed = fullBed.h

    print("Cylinder Radius (with rCylDelta):", R)
    print("Cylinder Height (full):", h)
    print("Bed Height (zmax - zmin):", hBed)

    # print("Cylinder Volume:", fullBed.CylinderVolume)
    # print("Packed Bed Volume:", fullBed.volume())
    # print("nBeads: ", len(fullBed.beads))

    cylvol = pi*R**2*h

    print("nBeads: ", len(fullBed.beads))
    print("Cylinder Volume:", cylvol)
    print("Packed Bed Volume:", fullBed.volume())
    print("Column Porosity:", 1-(fullBed.volume()/cylvol))
    print("Bed Porosity:", 1-(fullBed.volume()/fullBed.CylinderVolume))

    # # bridgeOffsetRatio = 0.95
    # bridgeOffsetRatio = sqrt(1 - relativeBridgeRadius**2)
    # bridgeTol = 0.04 * meshScalingFactor

    # addedBridgeVol, removedBridgeVol = bridgeVolumes(fullBed.beads, bridgeTol, relativeBridgeRadius, bridgeOffsetRatio)
    # print("Bridge Added Volume:", addedBridgeVol)
    # print("Capped Removed Volume:", removedBridgeVol)

    # sys.exit(0)

    print("rmin:", fullBed.rmin)
    print("rmax:", fullBed.rmax)
    print("ravg:", fullBed.ravg)

    volFrac, avgRads = histo([bead.r for bead in fullBed.beads], filename='psdtotal')

    _,BINS = np.histogram([bead.r for bead in fullBed.beads], bins=NBINS)

    print("\n--- Full Bed Histogram ---")
    print('vol_frac:\n', volFrac,)
    print('mean_radii:\n', avgRads)
    print("----\n")

    # pickler(volFrac, os.path.expanduser('~/fullbedvolfrac.pickle'))
    # pickler(avgRads, os.path.expanduser('~/fullbedavgRads.pickle'))

    nRegions = NREGIONS
    nShells = nRegions + 1 #Including r = 0
    rShells = []

    # shellType = 'EQUIVOLUME'
    shellType = 'EQUIDISTANT'

    if shellType == 'EQUIVOLUME':
        for n in range(nShells):
            rShells.append(R * sqrt(n/nRegions))
    elif shellType == 'EQUIDISTANT':
        for n in range(nShells):
            rShells.append(R * (n/nRegions))

    # print("rShells:", rShells)

    volRegions = [0] * nRegions

    ## Multiprocessing code.
    ##      Create a partial function of volShellRegion(beads, rShells, i) --> parfunc(i)
    ##      map each 'i' to each process
    pool = Pool()
    parfunc = partial(volShellRegion, fullBed.beads, rShells)
    # volRegions = pool.map(parfunc, range(nRegions))
    volRegions, radsRegion = zip(*pool.map(parfunc, range(nRegions)))
    pool.close()
    pool.join()

    # print(volShellRegion(fullBed.beads, rShells, 0))
    volRegions = [float(item) for item in volRegions]

    volCylRegions = [pi * hBed * (rShells[i+1]**2 - rShells[i]**2) for i in range(nRegions)]
    porosities = [ float(1-n/m) for n,m in zip(volRegions, volCylRegions) ]
    avg_shell_radii = [ (rShells[i] + rShells[i+1])/2 for i in range(nRegions) ]

    # for i in range(nRegions):
    #     print(avg_radius[i], volRegions[i], volCylRegions[i], porosities[i])

    print("\n--- Radial Porosity Distribution in Bed ---")
    # print("avg shell radii:\n", avg_shell_radii)
    print("porosities:\n", porosities)
    print("---\n")

    ######## pickler(avg_shell_radii, '~/avgshellradii.pickle')
    # pickler(porosities, os.path.expanduser('~/porosities.pickle'))

    volFracs = []
    for rads in radsRegion:
        volFrac, avgRads = histo([float(x) for x in rads], bins=BINS)
        volFracs.extend(volFrac)

    ## flatten volFracs
    # volFracs = [ v for sublist in volFracs for v in sublist ]

    print("\n--- Particle Distribution per Shell ---")
    print("vol_frac length (NRAD * NPARTYPE) =", len(volFracs))
    print("avg_radii length (NPARTYPE) =", len(avgRads))
    print("vol_frac =", volFracs)
    print("avg_radii = ", avgRads)
    print("---\n")

    # pickler(volFracs, os.path.expanduser( '~/volfracs.pickle') )
    # pickler(avgRads, os.path.expanduser( '~/avgrads.pickle') )

    plotter(avg_shell_radii, porosities, '', 'por_rad.pdf')
    csvWriter('por_rad.csv', avg_shell_radii, porosities)


def volShellRegion(beads, rShells, i):
    """
    Find the intersection volumes between rShells[i] & rShells[i+1]
    """
    volShell=0
    radsShell=[]
    for bead in beads:
        volBead = volBeadSlice(bead, rShells[i], rShells[i+1])
        volShell = volShell + volBead
        radBead = pow(volBead/(4.0/3.0*pi), 1.0/3.0)
        if radBead != 0:
            radsShell.append(radBead)
    return volShell, radsShell

# def radsShellRegion(beads, rShells, i):
#     """
#     In order to calculate histogram for beads in individual shells.
#     > 1. Calculate Volumes of each bead in given shell,
#     > 2. Extrapolate "radius" from each volume, even sliced ones.
#     > 3. Return list of radii to be used by histo()
#     """
#     radsShell=[]
#     for bead in beads:
#         volBead = volBeadSlice(bead, rShells[i], rShells[i+1])
#         radBead = pow(volBead/(4.0/3.0*pi), 1.0/3.0)
#         radsShell.append(radBead)
#     return radsShell

def bridgeVolumes(beads, bridgeTol, relativeBridgeRadius, bridgeOffsetRatio):
    """
    Find the total volume of the bridges between beads
    """
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
    # print(sys.version)
    # print(__doc__)
    main()
