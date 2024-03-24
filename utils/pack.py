#!/usr/bin/env python3

"""
### POTENTIALLY DEPRECATED: SEE pymesh/pack-info.py for meshes generated with pymesh instead!
## See WARNING below

@name: pack.py
@desc: Read binary packing.xyzd files, and generate porosity profiles.
@theory: Cylinder/Sphere intersection volume equations from http://dx.doi.org/10.1016/s1385-7258(61)50049-2
@usage: ./pack.py <packing.xyzd> <zBot after scaling> <zTop after scale> <scaling factor>
        zBot: bottom limit of slice to look for beads by center point.
        zTop: top limit of slice ...
        scaling factor: two packing.xyzd might not have the same dimensions. This helps fix that.

@NOTE: Doesn't work with porosity-controlled genmesh output meshes yet!!.
@NOTE: porosity is calculated for the cylinder length provided, not for void spaces yet
@NOTE: In this script, we only deal with radial variation, so other multiplexes (see cadet 2DGRM doc) are ignored
"""

# TODO: Handle rho == eta edge cases
# TODO: Auto handle scaling_factor: (Updatebounds, scale to fit Cyl Radius = 5)
# TODO: Better parallelization?

# TODO: Write output plots as pdf: 1. Histogram, 2. Radial porosity profile

import sys
import struct
import itertools
import numpy as np
import argparse
from matplotlib import pyplot as plt
from math import asin,sqrt,pi
from mpmath import ellipk, ellipe, ellipf, nstr
from multiprocessing import Pool
from functools import partial
import pickle
import csv

# NPARTYPE    = 10                   ## NPARTYPE in CADET == number of bins to sort particles by size
# NRAD = 20                          ## NRAD in CADET == Number of radial shells

def csvWriter(filename, x, y):
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

    def surface_area(self):
        return 4 * np.pi * self.r**2

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

    def surface_area(self):
        return sum([bead.surface_area() for bead in self.beads])

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
        z = []

        for bead in self.beads:
            xpr.append(bead.x + bead.r)
            xmr.append(bead.x - bead.r)
            ypr.append(bead.y + bead.r)
            ymr.append(bead.y - bead.r)
            zpr.append(bead.z + bead.r)
            zmr.append(bead.z - bead.r)
            z.append(bead.z)

        radList = [ bead.r for bead in self.beads ]
        self.rmax = max(radList)
        self.rmin = min(radList)
        self.ravg = sum(radList)/len(radList)

        self.bound_zbot = min( bead.z for bead in self.beads )

        self.xmax = max(xpr)
        self.ymax = max(ypr)
        self.ymin = min(ymr)
        self.xmin = min(xmr)
        self.zmax = max(zpr)
        self.zmin = min(zmr)

        self.dx = self.xmax - self.xmin
        self.dy = self.ymax - self.ymin
        self.dz = self.zmax - self.zmin

        self.R = max((self.xmax-self.xmin)/2, (self.ymax-self.ymin)/2) ## Similar to Genmesh
        self.h = self.zmax - self.zmin
        self.CylinderVolume = pi * self.R**2 * self.h

    def moveBedtoCenter(self):
        """
        Translate bed center to origin of coordinate system.
        """
        self.updateBounds()
        offsetx = -(self.xmax + self.xmin)/2
        offsety = -(self.ymax + self.ymin)/2
        offsetz = -(self.bound_zbot)
        for bead in self.beads:
            bead.x = bead.x + offsetx
            bead.y = bead.y + offsety
            bead.z = bead.z + offsetz
        self.updateBounds()

    def prune_to_volume(self, target_volume:float, eps:float = 1e-3): 
        """
        Prune packed bed of beads to reach a target volume
        """
        self.updateBounds()

        if self.volume() < target_volume: 
            print("Cannot prune packed bed! current volume < target volume")
            print(f"{self.volume() = }")
            print(f"{target_volume = }")
            sys.exit(-1)

        delta_volume = self.volume() - target_volume

        print(f"{self.volume() = }")
        print(f"{target_volume = }")
        print(f"{delta_volume = }")

        while delta_volume/target_volume > eps: 

            ## WARNING: This may not be the same logic used to create meshes
            # del_zone_beads = list(filter(lambda b: b.z < self.zmin + self.rmax, self.beads)) + list(filter(lambda b: b.z > self.zmax - self.rmax, self.beads)) 
            print("NOTE: Selecting beads in the end zone only for deletion.")
            del_zone_beads = list(filter(lambda b: b.z > self.zmax - self.rmax, self.beads)) 
            out = min(del_zone_beads, key=lambda b: abs(b.volume() - delta_volume))

            self.beads.remove(out)

            print(f"Deleting bead {out} with volume = {out.volume()}")

            self.updateBounds()
            delta_volume = self.volume() - target_volume

        print(f"{self.volume() = }")
        print(f"{target_volume = }")
        print(f"{delta_volume = }")
        print(f"{len(self.beads) = }")
        self.updateBounds()
        self.print_bounds()

    def print_bounds(self): 
        dic = {
                'xmin': self.xmin,
                'xmax': self.xmax,
                'ymin': self.ymin,
                'ymax': self.ymax,
                'zmin': self.zmin,
                'zmax': self.zmax,
                'rmin': self.rmin,
                'rmax': self.rmax,
                'ravg': self.ravg,
                'R': self.R,
                'h': self.h,
                'volume': self.volume()
                }
        print(dic)


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
    bins = kwargs.get('bins', 1)

    V=[4*np.pi*x*x*x/3 for x in radii]

    ## Dump into bins by weight of each bead's volume
    ## h (height of histogram bar) is then a representation of
    ## the volume of beads present at a certain radius (partype).
    ## if density==true weights are normalized
    h,e = np.histogram(radii, bins=bins, density=True, weights=V)


    ## Find the volume fraction at each point
    frac=[x/sum(h) for x in h]
    # print(sum(frac))

    ## Find means of each bin from the edges (e)
    w=2
    avg=np.convolve(e, np.ones(w), 'valid') / w

    if filename:
        with plt.style.context(['science']):
            # matplotlib.rcParams['font.sans-serif'] = "Verdana"
            # matplotlib.rcParams['font.family'] = "sans-serif"

            fig, ax = plt.subplots()
            ax.hist(radii, bins=bins)

            # ax.set(title=filename)
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
    ap.add_argument("-z", "--zlimits", type=float, nargs=2, help="bottom/top z value for packing slice. Limits based on bead centers.")
    ap.add_argument("-psf", "--pre-scaling-factor", type=float, default=1, help="Scaling factor to account for differences in mono/poly packing data")
    ap.add_argument("-msf", "--mesh-scaling-factor", type=float, default=1e-4, help="Post meshing scaling factor")
    ap.add_argument("-rf", "--r-factor", type=float, default=1, help="Bead radius shrinking factor")

    ap.add_argument("-c", "--container", default='cylinder', choices=['cylinder', 'box'], help="type of container")
    ap.add_argument("-ccsa", "--container-cross-section-area", type=float, help="Container cross section area")

    ap.add_argument("--nrad", type=int, default=1, help="NRAD, number of radial shells in 2D model")
    ap.add_argument("--npartype", type=int, default=1, help="NPARTYPE, number of bins to sort particles by size")
    ap.add_argument("-st"  , "--shelltype", choices = ['EQUIDISTANT', 'EQUIVOLUME'], default='EQUIDISTANT', help="Shell discretization type")
    ap.add_argument("--prune-to-volume", type=float, default=None, help="Shell discretization type")
    ap.add_argument("--eps", type=float, default=1e-3, help="Epsilon for prune-to-volume check")

    ap.add_argument("-v", "--vartype", help="type of variable stored (d | f | i) in packing file", default='d')
    ap.add_argument("-e", "--endianness", help="> or <", default='<')

    ap.add_argument("-d", "--dry-run", help="Do not calculate vol_frac, porosities and mean_radii", action='store_true')
    ap.add_argument("-o", "--output-prefix", help="prefix to output for radial porosity plot and PSD histogram", required=True)

    ap.add_argument("--inlet", type=float, default=2.5, help="Void space before zlimits.")
    ap.add_argument("--outlet", type=float, default=2.5, help="Void space after zlimits.")

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

        endianness = '<'
        packingPrecision = int(infiledict['packingPrecision'])
        if packingPrecision == 8:
            vartype = 'd'
        else:
            vartype = 'f'
    else:
        packing = args['packing']
        vartype = args['vartype']
        endianness = args['endianness']
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
    inlet = args['inlet'] * meshScalingFactor
    outlet = args['outlet'] * meshScalingFactor

    fullBed = PackedBed()

    # dataformat = "<f" ## For old packings with little endian floating point data. Use <d for new ones
    dataformat = endianness + vartype
    arr = bin_to_arr(packing, dataformat)
    for chunk in grouper(arr,4):
        if (chunk[2] >= zBot/scaling_factor) and (chunk[2] <= zTop/scaling_factor):
            x = chunk[0] * scaling_factor * meshScalingFactor
            y = chunk[1] * scaling_factor * meshScalingFactor
            z = chunk[2] * scaling_factor * meshScalingFactor
            r = chunk[3]/2 * scaling_factor * rFactor * meshScalingFactor
            fullBed.add(Bead(x, y, z, r))

    fullBed.moveBedtoCenter()

    if(args['prune_to_volume']):
        fullBed.prune_to_volume(args['prune_to_volume'], eps=args['eps'])

    ## WARNING: Note the way the column length and radius are calculated,
    ## As in genmesh, we set a zTop and zBot, and add inlet + outlet to the ends
    ## In pymesh, however, we define container dimensions directly.
    ## pymesh/pack-info.py should be more apt for meshes generated via pymesh.
    R = fullBed.R + rCylDelta ## Adding Rcyldelta
    h = (zTop - zBot)*meshScalingFactor  + inlet + outlet
    hBed = fullBed.h

    print("xmax = {}, xmin = {}".format(fullBed.xmax, fullBed.xmin))
    print("ymax = {}, ymin = {}".format(fullBed.ymax, fullBed.ymin))
    print("zmax = {}, zmin = {}".format(fullBed.zmax, fullBed.zmin))

    print("Cylinder Radius (with rCylDelta):", R)
    print("Container Height (full):", h)
    print("Bed Height (zmax - zmin):", hBed)
    print("Inlet length till bead zBot: ", inlet)
    print("Outlet length till bead zBot: ", outlet)
    print("Bed surface area: ", fullBed.surface_area())

    # print("Cylinder Volume:", fullBed.CylinderVolume)
    # print("Packed Bed Volume:", fullBed.volume())
    # print("nBeads: ", len(fullBed.beads))


    if args['container'] == 'cylinder': 
        container_cross_section_area = pi*R**2
    elif args['container'] == 'box': 
        container_cross_section_area = fullBed.dx * fullBed.dy

    if args['container_cross_section_area']: 
        container_cross_section_area = args['container_cross_section_area'] * args['mesh_scaling_factor']**2

    container_volume = container_cross_section_area * h
    container_bedlength_volume = container_cross_section_area * hBed

    print("nBeads: ", len(fullBed.beads))
    print("Container Volume:", container_volume)
    print("Packed Bed Volume:", fullBed.volume())
    print("Container Bed-Length Volume:", container_bedlength_volume)
    print("Column Porosity:", 1-(fullBed.volume()/container_volume))
    print("Bed Porosity:", 1-(fullBed.volume()/container_bedlength_volume))


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

    if args['dry_run']:
        sys.exit(0)

    with open(args['output_prefix'] + '_beads.csv', 'w') as f:
        writer = csv.writer(f)
        writer.writerows([[bead.r] for bead in fullBed.beads])

    volFrac, bin_radii = histo([bead.r for bead in fullBed.beads], filename=args['output_prefix'] + '_psd', bins=args['npartype'])

    # If bins is an int, it defines the number of equal-width bins in the given range
    # (10, by default). If bins is a sequence, it defines a monotonically increasing
    # array of bin edges, including the rightmost edge, allowing for non-uniform bin widths.

    # Generate BINS, based on all the beads in the whole column. This will be used later
    # to categorize beads within shells into the same BINS
    _,BINS = np.histogram([bead.r for bead in fullBed.beads], bins=args['npartype'])

    print("\n--- Full Bed Histogram ---")
    print('vol_frac:\n', volFrac,)
    print('mean_radii:\n', bin_radii)
    print("----\n")

    nRegions = args['nrad']
    nShells = nRegions + 1 #Including r = 0
    rShells = []

    shellType = args['shelltype']

    if shellType == 'EQUIVOLUME':
        for n in range(nShells):
            rShells.append(R * sqrt(n/nRegions))
    elif shellType == 'EQUIDISTANT':
        for n in range(nShells):
            rShells.append(R * (n/nRegions))

    # print("rShells:", rShells)

    total_beads_volume_per_shell = [0] * nRegions

    ## Multiprocessing code.
    ##      Create a partial function of volShellRegion(beads, rShells, i) --> parfunc(i)
    ##      map each 'i' to each process
    pool = Pool()
    parfunc = partial(volShellRegion, fullBed.beads, rShells)
    # volRegions = pool.map(parfunc, range(nRegions))
    total_beads_volume_per_shell, radii_beads_per_shell = zip(*pool.map(parfunc, range(nRegions)))
    pool.close()
    pool.join()

    # print(volShellRegion(fullBed.beads, rShells, 0))

    total_beads_volume_per_shell = np.array(total_beads_volume_per_shell).astype(np.float64)
    radii_beads_per_shell = [ np.array(item).astype(np.float64) for item in radii_beads_per_shell ]

    volCylRegions_bed = [pi * hBed * (rShells[i+1]**2 - rShells[i]**2) for i in range(nRegions)]
    volCylRegions_column = [pi * h * (rShells[i+1]**2 - rShells[i]**2) for i in range(nRegions)]

    porosities_bed = [ float(1-n/m) for n,m in zip(total_beads_volume_per_shell, volCylRegions_bed) ]
    porosities_column = [ float(1-n/m) for n,m in zip(total_beads_volume_per_shell, volCylRegions_column) ]

    avg_shell_radii = [ (rShells[i] + rShells[i+1])/2 for i in range(nRegions) ]

    print("\n--- Radial Porosity Distribution in Bed ---")
    print("col_porosity_bed:\n", porosities_bed)
    print("---\n")

    print("\n--- Radial Porosity Distribution in FULL COLUMN ---")
    print("col_porosity:\n", porosities_column)
    print("---\n")

    ## Get histogram data: volume fractions and radii, for each shell
    ## bin_radii is the list of mean bin radii for each shell, which is set to BINS
    volFracs = []
    for rads in radii_beads_per_shell:
        volFrac, bin_radii = histo([float(x) for x in rads], bins=BINS)
        volFracs.extend(volFrac)

    print("\n--- Particle Distribution per Shell ---")
    print("vol_frac length (NRAD * NPARTYPE) =", len(volFracs))
    print("avg_radii length (NPARTYPE) =", len(bin_radii))
    print("par_type_volfrac =", volFracs)
    print("par_radius = ", bin_radii)
    print("---\n")

    plotter(avg_shell_radii, porosities_bed, '', args['output_prefix'] + '_bedpor_rad.pdf')
    plotter(avg_shell_radii, porosities_column, '', args['output_prefix'] + '_colpor_rad.pdf')
    csvWriter(args['output_prefix'] + '_bedpor_rad.csv', avg_shell_radii, porosities_bed)
    csvWriter(args['output_prefix'] + '_colpor_rad.csv', avg_shell_radii, porosities_column)

def volShellRegion(beads, rShells, i):
    """
    Find the intersection volumes between rShells[i] & rShells[i+1]

    @input: beads, shell_radii, index of shell
    @output:
        - total volume of all particles within the i'th shell
        - list of all radii based on intersected volumes.

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
