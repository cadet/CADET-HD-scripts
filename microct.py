#!/usr/bin/env python3

"""
MicroCT image analysis on collections of images
@task: Average over multiple slices
@task: is it possible that the variation is experimental systematic error due to light dissipating in the medium?
@task: xns brinkman simulations with radial porosity profile
"""


import matplotlib.pyplot as plt
from skimage.filters import threshold_mean, threshold_yen, threshold_otsu
from skimage.morphology import binary_dilation, binary_erosion
from skimage import io
import numpy as np
import argparse
import csv

def plot_image_threshold(img, threshold, binary, filename=None):
    fig, axes = plt.subplots(ncols=3, figsize=(8, 2.5))
    ax = axes.ravel()
    ax[0] = plt.subplot(1, 3, 1, adjustable='box')
    ax[1] = plt.subplot(1, 3, 2)
    ax[2] = plt.subplot(1, 3, 3, sharex=ax[0], sharey=ax[0], adjustable='box')

    ax[0].imshow(img, cmap=plt.cm.gray) # type: ignore
    ax[0].set_title('Original')
    ax[0].axis('off')

    ax[1].hist(img.ravel(), bins=256)
    ax[1].set_title('Histogram')
    ax[1].axvline(threshold, color='r')

    ax[2].imshow(binary, cmap=plt.cm.gray) # type: ignore
    ax[2].set_title('Thresholded')
    ax[2].axis('off')

    if filename:
        fig.savefig(filename)

    plt.show()

def plot_image(img, filename=None):
    plt.imshow(img, cmap=plt.cm.gray) # type: ignore

    if filename:
        plt.savefig(filename)
    plt.show()

def die():
    import sys
    sys.exit()

def driver(img, half_width, threshold='mean',  filename='output'):
    """
    Calculate porosity profile for a single image
    """

    if threshold == 'mean':
        threshold = threshold_mean(img)
    elif threshold == 'otsu':
        threshold = threshold_otsu(img)
    elif threshold == 'yen':
        threshold = threshold_yen(img)
    else:
        threshold = float(threshold)

    print(f"{threshold=}")
    binary = img > threshold

    ## arr is used for further operations. If set to img, we can use the threshold
    ## calculated directly on it in later operations.
    ## If we use binary, we have to set the threshold to 0 or 1 or something
    arr = np.array(img)
    # arr = binary

    ## Calculate xy bounds of the image
    x = np.arange(0, arr.shape[1])
    y = np.arange(0, arr.shape[0])
    xmax, ymax = arr.shape
    print(arr.shape)
    # print(arr[0][0])

    dilation1 = binary_dilation(binary)
    erosion1 = binary_erosion(dilation1)
    erosion2 = binary_erosion(erosion1)
    dilation2 = binary_dilation(erosion2)

    plot_image(binary, filename=filename + '_binary.pdf')
    plot_image(dilation2, filename=filename + '_cleaned.pdf')

    ## Find limits and center of the column.
    yw,xw = np.where(arr >= threshold) ## Find bright pixels -> beads

    xmin = min(xw)
    xmax = max(xw)
    ymin = min(yw)
    ymax = max(yw)

    ## Print column limits
    print("Column X: [{}, {}]".format(xmin, xmax))
    print("Column Y: [{}, {}]".format(ymin, ymax))

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
    #arr[mask] = 3
    ## arr[mask] = 1
    #colorplot(x,y,arr,"microct-sample-ring.png")
    #sys.exit()

    ## Output lists
    out_rad = []
    out_por = []

    out_rad_por_tups = []

    zarr = np.zeros((arr.shape[0], arr.shape[1]))

    for r in range(half_width,rad-half_width,2*half_width):
        # print(r-half_width, r, r+half_width)

        ## Find all pixels within ring of radius: r-half_width to r+half_width
        mask = np.logical_and( (x[np.newaxis,:]-xc)**2 + (y[:,np.newaxis]-yc)**2 <= (r+half_width)**2, (x[np.newaxis,:]-xc)**2 + (y[:,np.newaxis]-yc)**2 >= (r-half_width)**2)

        total = arr[mask].size ## Count number of pixels in given ring
        pores = np.where(arr[mask] < threshold)[0].size ## count number of dark pixels -> voids
        # print(pores, total)

        # ## To check if I'm double sampling any pixels
        # ## Turns out, yes, but very few to actually matter
        # zarr[mask] = zarr[mask] + 1

        # out_rad.append(r)
        # out_por.append(pores/total)
        out_rad_por_tups.append((r,pores/total))

    out_por_avg = np.average(out_por)

    eps = 1e-3
    out_rad_por_tups = [ tup for tup in out_rad_por_tups if tup[0] > 50 and tup[1] < 1-eps ]
    por_avg = sum(map(lambda x: x[1], out_rad_por_tups)) / len(out_rad_por_tups)

    ## Plotting
    # fig = plt.figure()
    # plt.plot(out_rad, out_por,linewidth=2)
    # plt.xlabel("Radial distance")
    # plt.ylabel("Porosity")
    # plt.show()
    # fig.savefig("image_rad_por.png")

    ## Full disk mask
    # fulldiskmask = ((x[np.newaxis,:]-xc)**2 + (y[:,np.newaxis]-yc)**2 <= (rad)**2 )
    # zarr[fulldiskmask] = 1
    # colorplot(x,y,zarr, "disk.png")
    # total = arr[fulldiskmask].size ## Count number of pixels in given ring
    # pores = np.where(arr[fulldiskmask] < threshold)[0].size ## count number of dark pixels -> voids
    # print("Average porosity of current slice: {avg}".format(avg=pores/total))

    with open (filename + '_rad_por_tups.csv', 'w') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerows(out_rad_por_tups)

    with open (filename + '_rad_por_avg.csv', 'w') as f:
        writer = csv.writer(f, delimiter=',')

        ## Average including "bad" areas
        # writer.writerow((out_rad[0], out_por_avg))
        # writer.writerow((out_rad[-1], out_por_avg))

        ## Average after snipping "bad" areas
        writer.writerow((0, por_avg))
        writer.writerow((rad, por_avg))

def main():

    ap = argparse.ArgumentParser()
    ap.add_argument('-w', '--width', type=int, default=1, help='Half width of the sampling ring in pixels.')
    ap.add_argument('-t', '--threshold', default='mean', help='Thresholding function or value')
    # ap.add_argument('-a', '--average', action='store_true', help="Create an average plot over multiple files.")
    ap.add_argument('FILES', nargs='*', help='Image files to process.')
    args = vars(ap.parse_args())

    half_width    = args['width']     ## Half ring-width to count pixels at a given radius

    coll = io.collection.ImageCollection(args['FILES'])
    arr = coll.concatenate()
    print(arr.shape)
    num_files = arr.shape[0]

    ## NOTE: Can parallelize
    for i in range(num_files):
        driver( arr[i,:,:], args['width'], args['threshold'], filename='output_{:04d}'.format(i) )

if __name__ == "__main__":
    main()
