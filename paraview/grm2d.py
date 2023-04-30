#!/usr/bin/env pvpython

from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()
from vtkmodules.numpy_interface import dataset_adapter as dsa
import vtk.util.numpy_support as ns #type:ignore
from math import sqrt
import numpy as np
import struct

import argparse
import os
import csv

def main():


    ap = argparse.ArgumentParser()

    ap.add_argument("--ncol", default=1, type=int, help="Number of axial sections")
    ap.add_argument("--nrad", default=1, type=int, help="Number of radial sections")
    ap.add_argument("-st"  , "--shelltype", choices = ['EQUIDISTANT', 'EQUIVOLUME'], default='EQUIDISTANT', help="Shell discretization type")

    ap.add_argument("-f", "--filetype", default='pvtu', choices=['xdmf', 'vtu', 'vtk', 'pvtu'], help="filetype: xdmf | vtu | vtk | pvtu")
    ap.add_argument("-s", "--scalars" , nargs='*' , help="Scalars to consider. (Previously colorvars).")
    ap.add_argument("FILES", nargs='*', help="files..")

    args = vars(ap.parse_args())

    if len(args['FILES']) == 0:
        filetype = args['filetype']

        try:
            args['FILES'] = sorted([ file for file in os.listdir(".") if file.endswith(filetype) ], key=lambda x: int(x.split('.')[0].split('_')[-1]))
        except:
            print("Not sorting files.")
            args['FILES'] = [ file for file in os.listdir(".") if file.endswith(filetype) ]

        if len(args['FILES']) == 0:
            print("Didn't find", filetype, "files in current folder.")
            sys.exit(-1)
    else:
        fileExtensions = set([os.path.splitext(infile)[1] for infile in args['FILES']])
        if len(fileExtensions) > 1:
            print("Mixed File Formats Given!")
            sys.exit(-1)
        filetype = fileExtensions.pop().replace('.', '')

    for key in args:
        print(key + ': ', args[key])

    reader=None
    if filetype == 'xdmf':
        reader = XDMFReader(FileNames=args['FILES'])
    elif filetype == 'vtu':
        reader = XMLUnstructuredGridReader(FileName=args['FILES'])
    elif filetype == 'pvtu':
        reader = XMLPartitionedUnstructuredGridReader(FileName=args['FILES'])
    elif filetype == 'vtk':
        reader = LegacyVTKReader(FileNames=args['FILES'])
    else:
        print("Unsupported File Format!")
        raise(ValueError)

    timeKeeper = GetTimeKeeper()

    timeArray = reader.TimestepValues
    nts = len(timeArray) or 1

    args['scalars'] = args['scalars'] or reader.PointArrayStatus

    scalars = args['scalars']
    shellType = args['shelltype']

    ## Split into axial columns
    ## Split into cylindrical columns
    ## Integrate

    print("N =", servermanager.ActiveConnection.GetNumberOfDataPartitions())
    print("n =", servermanager.vtkProcessModule.GetProcessModule().GetPartitionId())

    import sys; sys.exit()

    # view = GetActiveViewOrCreate('RenderView')
    # display = Show(object, view)
    # display.Representation = args['display_representation']
    # (xmin,xmax,ymin,ymax,zmin,zmax) = GetActiveSource().GetDataInformation().GetBounds()
    # radius = ((xmax-xmin) + (ymax-ymin)) / 4
    # length = zmax - zmin
    # print("Length: {}, Radius: {}".format(length, radius))

    # # nCol = 10              # Number of axial regions
    # # nRad = 5               # Number of radial regions
    # nCol = args['ncol']
    # nRad = args['nrad']

    # # dx = length/nCol
    # # dr = radius/nRad

    # colEdges = np.linspace(zmin, zmax, nCol+1)
    # radEdges = np.linspace(0,radius,nRad+1) if args['shelltype'] == 'EQUIDISTANT' else list(x/nRad * radius for x in range(nRad+1))


    # # Hide(object, view)
    # nColEdgeFractions = np.linspace(0,1,nCol+1)
    # nRadEdgeFractions = np.linspace(0,1,nRad+1)

    # ## Output vector. Should contain (nts X nCol X nRad X nScalar)
    # grm2d_output = []

    # ## Hack to remove previous file
    # arr_to_bin([], 'grm2doutput_unpacked.bin', 'd')

    # for timestep in range(nts):

    #     timeKeeper.Time = timestep
    #     # object.UpdatePipeline(reader.TimestepValues[timestep])
    #     # object.UpdatePipeline(timeArray[timestep])
    #     object.UpdatePipeline()

    #     grm2d_timestep_output = []

    #     print("--> TS: {}".format(timestep))

    #     # for leftEdge, rightEdge in zip(colEdges[:-1], colEdges[1:]):
    #     for leftEdge, rightEdge in zip(nColEdgeFractions[:-1], nColEdgeFractions[1:]):
    #         print("  |--> Col: {}/{}".format(np.where(nColEdgeFractions == rightEdge)[0][0],nCol))
    #         SetActiveSource(object)
    #         # print('[{}, {}]'.format(leftEdge, rightEdge))

    #         clipLeftArgs = { 'project' : ['clip', 'Plane', leftEdge , '-z'] }
    #         clipRightArgs = { 'project' : ['clip', 'Plane', rightEdge, '+z'] }

    #         clipLeft = project(object, clipLeftArgs)
    #         clipRight = project(clipLeft, clipRightArgs)

    #         radAvg = []

    #         for radIn, radOut in zip(radEdges[:-1], radEdges[1:]):
    #             radAvg.append( (radIn + radOut) / 2 )
    #             # print('--> [{}, {}]: {}'.format(radIn, radOut, (radIn+radOut)/2))

    #             print('    |--> Rad: {}/{}'.format(np.where(radEdges == radOut)[0][0], nRad))

    #             clipOuter = Clip(Input=clipRight)
    #             clipOuter.ClipType = 'Cylinder'
    #             clipOuter.ClipType.Axis = [0.0, 0.0, 1.0]
    #             clipOuter.ClipType.Radius = radOut
    #             Hide3DWidgets(proxy=clipOuter.ClipType)

    #             # renderView1 = GetActiveViewOrCreate('RenderView')
    #             # projectionDisplay = Show(clipOuter, renderView1)
    #             # projectionDisplay.Representation = 'Surface'
    #             # # projectionDisplay.Representation = 'Surface With Edges'
    #             # renderView1.OrientationAxesVisibility = int(axisVisible)
    #             # projectionDisplay.RescaleTransferFunctionToDataRange()

    #             clipInner = Clip(Input=clipOuter)
    #             clipInner.ClipType = 'Cylinder'
    #             clipInner.ClipType.Axis = [0.0, 0.0, 1.0]
    #             clipInner.ClipType.Radius = radIn
    #             clipInner.Invert = 0

    #             integrated_scalars = integrate(clipInner, args['scalars'], normalize='Volume')
    #             # print('---->', integrated_scalars[0])
    #             grm2d_timestep_output.extend(integrated_scalars[0])

    #             Delete(clipInner)
    #             Delete(clipOuter)

    #         Delete(clipLeft)
    #         Delete(clipRight)

    #     grm2d_output.extend(grm2d_timestep_output)
    #     # arr_to_bin_unpacked(grm2d_timestep_output, 'grm2doutput_unpacked.bin', 'd', mode='a')

    # # print(grm2d_output)
    # # arr_to_bin(grm2d_output, 'grm2doutput.bin', 'd')
    # print("DONE!")
    # arr_to_bin_unpacked(grm2d_output, 'grm2doutput_unpacked.bin', 'd')

if __name__ == "__main__":
    main()
