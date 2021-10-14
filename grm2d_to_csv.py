#!/usr/bin/env python3

import argparse
import struct
import numpy as np

def csvWriter(filename, data):
    import csv
    with open(filename, 'w') as f:
        writer = csv.writer(f)
        # writer.writerows(zip(x, y))
        writer.writerows(zip(data))

def dataformatsize(dataformat):
    vartype = dataformat[1]
    datasize = 0
    if vartype == 'd':
        datasize = 8
    elif vartype == 'f':
        datasize = 4
    elif vartype == 'i':
        datasize = 4
    return datasize


def bin_to_arr(filename, dataformat, skip=0, skiprows=0, nrows=0, ncols=0):
    datasize = dataformatsize(dataformat)
    arr = []

    with(open(filename, 'rb')) as input:
        input.seek(skip * nrows * ncols * datasize + skiprows * ncols * datasize, 0)
        myiter = struct.iter_unpack(dataformat, input.read())

        for i in myiter:
            arr.append(i[0])

    return arr

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("FILES", nargs='*', help="files")
    ap.add_argument("--nts", type=int, help="num. timesteps")
    ap.add_argument("--ncol", type=int, help="num. column segments")
    ap.add_argument("--nrad", type=int, help="num. radial segments")
    ap.add_argument("--nscalar", type=int, default=1, help="num. scalars")
    ap.add_argument("-o", "--output", default="output.csv", help="Output file")
    ap.add_argument("-t", "--timesnap", nargs=3, type=int, help="[icol, irad, iscalar] time snapshot")
    ap.add_argument("-c", "--colsnap", nargs=3, type=int, help="[its, irad, iscalar] column length snapshot")
    ap.add_argument("-r", "--radsnap", nargs=3, type=int, help="[its, icol, iscalar] radial snapshot")
    ap.add_argument("--timeloop", action='store_true', help="write out conc_time data for each col_rad_scalar")
    ap.add_argument("--colloop", action='store_true', help="write out conc_col data for each time_rad_scalar")
    ap.add_argument("--radloop", action='store_true', help="write out conc_rad data for each time_col_scalar")
    args = vars(ap.parse_args())

    data = []

    for filename in args['FILES']:
        data.extend(bin_to_arr(filename, '<d'))

    nts = args['nts']
    ncol = args['ncol']
    nrad = args['nrad']
    nscalar = args['nscalar']

    ## nts, ncol, nrad, nscalar
    data = np.reshape(data, (nts, ncol, nrad, nscalar))

    ## NOTE: Saves to CSV conc_time data for each shell

    if args['timeloop']:
        for icol in range(ncol):
            for irad in range(nrad):
                for iscalar in range(nscalar):
                    csvWriter("timeloop_{}_{}_{}.csv".format(icol, irad, iscalar),data[:,icol,irad,iscalar])
    elif args['colloop']:
        for its in range(nts):
            for irad in range(nrad):
                for iscalar in range(nscalar):
                    csvWriter("colloop_{}_{}_{}.csv".format(its, irad, iscalar),data[its,:,irad,iscalar])
    elif args['radloop']:
        for its in range(nts):
            for icol in range(ncol):
                for iscalar in range(nscalar):
                    csvWriter("radloop_{}_{}_{}.csv".format(its, icol, iscalar),data[its,icol,:,iscalar])


    if args['timesnap']:
        icol = args['timesnap'][0]
        irad = args['timesnap'][1]
        iscalar = args['timesnap'][2]
        csvWriter(args['output'], data[:, icol, irad, iscalar])
    elif args['colsnap']:
        its= args['colsnap'][0]
        irad = args['colsnap'][1]
        iscalar = args['colsnap'][2]
        csvWriter(args['output'], data[its, :, irad, iscalar])
    elif args['radsnap']:
        its= args['radsnap'][0]
        icol= args['radsnap'][1]
        iscalar = args['radsnap'][2]
        csvWriter(args['output'], data[its, icol, :, iscalar])



if __name__ == "__main__":
    main()
