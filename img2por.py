#!/usr/bin/python3

"""
Basic image analysis script to calculate radial porosity profile from a radial snapshot, part of a micro-CT imaging.

Assumes beads => bright, voids = dark

Algorithm:
    - Enhance image contrast by enhanceFactor
    - Find column center using provided threshold
    - For every ring (radius=r from 0 to R, width=2*half_width),
        - calculate porosity as num(dark)/num(total) pixels

@usage: ./img2por.py <image.jpeg> -t <threshold> -w <ring-half-width> -e <enhanceFactor>
@note: for some reason, tiff images don't work with PIL that well, use JPEG.
"""

## TODO: Ensure that rings don't overlap.
## TODO: General refactoring and cleanup
## TODO: cmdline arg/flag to dump output images
## TODO: Debug mode to dump threshold boundary and ring width
## TODO: Implement axial averaging (in batches based on 2DGRM discretization)

from PIL import Image
from PIL import ImageEnhance
from PIL import ImageFilter
import numpy as np
import matplotlib.pyplot as plt
import argparse
import csv
import sys

# enhanceFactor = 1.5 ## Factor by which to enhance image contrast
# threshold     = 180 ## Brightness threshold to consider as beads/pores
# half_width    = 1   ## Half ring-width to count pixels at a given radius

ap = argparse.ArgumentParser()
ap.add_argument("-e", "--enhance", type=float, default=1.0, help="Enhance image contrast by this factor.")
ap.add_argument("-t", "--threshold", type=int, choices=range(256), required=True, help="Threshold to consider beads (bright) or void (dark).")
ap.add_argument("-w", "--width", type=int, default=1, help="Half width of the sampling ring in pixels.")
ap.add_argument("FILE", help="Image file to process.")
args = vars(ap.parse_args())

enhanceFactor = args['enhance']   ## Factor by which to enhance image contrast
threshold     = args['threshold'] ## Brightness threshold to consider as beads/pores
half_width    = args['width']     ## Half ring-width to count pixels at a given radius

im = Image.open(args['FILE'])
print(im.format, im.size, im.mode)

arr = np.array(im)
# print(arr.shape)

enh = ImageEnhance.Contrast(im)
ime = enh.enhance(enhanceFactor) # enhance contrast
# ime.show("ENHANCED")
# sys.exit()

# imb = ime.filter(ImageFilter.BLUR)
# imb2 = ime.filter(ImageFilter.MinFilter(3))
# imb3 = ime.filter(ImageFilter.MinFilter)
# imfe = ime.filter()
# imfe.show('blur', 'eog')
# sys.exit()

arr = np.array(ime)
x = np.arange(0, arr.shape[1])
y = np.arange(0, arr.shape[0])
xmax, ymax = arr.shape
# print(arr.shape)
# print(arr[0][0])


## Find limits and center of the column.
yw,xw = np.where(arr >= threshold) ## Find bright pixels -> beads
# print(xw.shape, yw.shape)

xmin = min(xw)
xmax = max(xw)
ymin = min(yw)
ymax = max(yw)

print(xmin, xmax)
print(ymin, ymax)

xrad = (xmax - xmin)//2
yrad = (ymax - ymin)//2
rad = (xrad + yrad) //2
print("Rads:", xrad, yrad, rad)

xc = (xmax + xmin) //2
yc = (ymax + ymin) //2
print("Centers:", xc, yc)


##Visual debug
## arr = np.zeros((arr.shape[0], arr.shape[1]))
#r = rad
#mask = np.logical_and( (x[np.newaxis,:]-xc)**2 + (y[:,np.newaxis]-yc)**2 <= (r+half_width)**2, (x[np.newaxis,:]-xc)**2 + (y[:,np.newaxis]-yc)**2 >= (r-half_width)**2)
## arr[mask] = 255
#arr[mask] = 0
#fig = plt.figure()
#plt.pcolormesh(x, y, arr)
#plt.colorbar()
#plt.show()
#fig.savefig("microct-image-400_sample_ring.png")

## Output lists
out_rad = []
out_por = []

for r in range(half_width,rad-half_width,2*half_width):
    # print(r-half_width, r, r+half_width)

    ## Find all pixels within ring of radius: r-half_width to r+half_width
    mask = np.logical_and( (x[np.newaxis,:]-xc)**2 + (y[:,np.newaxis]-yc)**2 <= (r+half_width)**2, (x[np.newaxis,:]-xc)**2 + (y[:,np.newaxis]-yc)**2 >= (r-half_width)**2)

    total = arr[mask].size ## Count number of pixels in given ring
    pores = np.where(arr[mask] <= threshold)[0].size ## count number of dark pixels -> voids

    out_rad.append(r)
    out_por.append(pores/total)

fig = plt.figure()
plt.plot(out_rad, out_por)
plt.show()
# fig.savefig("image_rad_por.png")

with open(args['FILE']+'_rad_por.csv', 'w') as f:
    writer = csv.writer(f, delimiter=',')
    writer.writerows(zip(out_rad, out_por))

