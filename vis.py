#!/bin/env -S pvpython --force-offscreen-rendering

"""
A Paraview script to aid the postprocessing of data. I wrote it mainly to generate images out of output files.

Usage by examples:
    > ./vis.py -h
    > pvpython --force-offscreen-rendering /path/to/vis.py **/*pvtu -s -c scalar_2 scalar_3 -p slice  ## for listed pvtu files, output scalar_2 and scalar_3 slice data
    > mpirun -np 4 pvbatch /path/to/vis.py -f pvtu -s                   ## looks for pvtu files in current folder, clips them, snapshots all scalars,
    > mpirun -np 4 pvbatch /path/to/vis.py **/*pvtu -d out-distributed  ## for listed pvtu files, apply d3 filter and write to out-distributed
    > mpirun -np 4 pvbatch /path/to/vis.py example-case_{150..300}.pvtu -s  ## Uses bash's sequence expansion

Usage notes:
    > don't forget to use a processing flag like -s or -d or -w
    > -c flag takes a list as input, make sure to use it AFTER filepaths
    > -f filetype is optional. Default is pvtu.
    > The program ALWAYS starts animation outputs with suffix 0000, regardless of the timestep data you've supplied. So be wary of using it twice in the same dir.
    > Assumes that timestep suffix is "_%d" to autodetect files in current directory

"""

## DONE: automate reader type
## DONE: distribute to parallel
## DONE: Fix flow animate problem
## DONE: Allow screenshotting ALL scalars in one run
## DONE: Parametrize view size geometry
## TODO: config.json file
## TODO: Parametrize background color
## TODO: fix timestep data (probably in mixd2pvtu)
## TODO: rotated views [Good for image resolution in long columns]
## TODO: After updatescalarbars() they are always visible even with -nsb

import argparse
from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()
from vtkmodules.numpy_interface import dataset_adapter as dsa
import vtk.util.numpy_support as ns
import numpy as np
from matplotlib import pyplot as plt
import os
import sys

def MIXDReader():
    pass

def mass_defect(reader, **kwargs):

    renderView1 = GetActiveViewOrCreate('RenderView')
    display = Show(reader, renderView1)

    (xmin,xmax,ymin,ymax,zmin,zmax) = GetActiveSource().GetDataInformation().GetBounds()

    zarr = []

    for zpos in np.linspace(1.0*zmin, zmax, 100):
        projection = Slice(Input=reader)
        Hide3DWidgets(proxy=projection.SliceType)

        projection.SliceType = 'Plane'
        projection.HyperTreeGridSlicer = 'Plane'
        projection.SliceOffsetValues = [0.0]

        projection.SliceType.Origin = [0.0, 0.0, zpos]
        projection.SliceType.Normal = [0.0, 0.0, 1.0]

        cellSize = CellSize(Input=projection)
        cellSize.ComputeSum = 1

        ## NOTE: POINT DATA TO CELL DATA BEFORE INTEGRATE?
        ## TODO: Get Bounding Box
        integrated = IntegrateVariables(Input=projection)

        intdata = servermanager.Fetch(integrated)
        intdata = dsa.WrapDataObject(intdata)
        value = intdata.PointData['scalar_3']
        zarr.append(value)

    zarr=[ns.vtk_to_numpy(x)[0] for x in zarr[2:]]
    print(zarr)
    plt.figure()
    plt.plot(zarr)
    plt.savefig('plot.pdf')
    # plt.show()


def snapshot(reader, **kwargs):
    projectionType = kwargs.get('projectionType', 'clip')
    colorVars = kwargs.get('colorVars', reader.PointArrayStatus) or reader.PointArrayStatus
    scalarBarVisible = kwargs.get('scalarBarVisible', True)
    axisVisible = kwargs.get('axisVisible', True)
    geometry = kwargs.get('geometry')

    animationScene1 = GetAnimationScene()
    timeKeeper1 = GetTimeKeeper()
    animationScene1.UpdateAnimationUsingDataTimeSteps()

    # try:
    #     ## Use last timestep as reference for creating color map
    #     animationScene1.AnimationTime = reader.TimestepValues[-1]
    #     timeKeeper1.Time = reader.TimestepValues[-1]
    # except:
    #     ## for files without time data
    #     animationScene1.AnimationTime = 0
    #     animationScene1.StartTime = 0
    #     animationScene1.EndTime = 0
    #     timeKeeper1.Time = 0

    projection = None

    if projectionType == 'clip':
        projection = Clip(Input=reader)
        Hide3DWidgets(proxy=projection.ClipType)
    elif projectionType == 'slice':
        projection = Slice(Input=reader)
        Hide3DWidgets(proxy=projection.SliceType)
    else:
        print("Invalid Projection Type!")
        sys.exit(-1)

    renderView1 = GetActiveViewOrCreate('RenderView')
    projectionDisplay = Show(projection, renderView1)
    projectionDisplay.Representation = 'Surface'
    projectionDisplay.SetScalarBarVisibility(renderView1, scalarBarVisible)
    renderView1.OrientationAxesVisibility = int(axisVisible)
    projectionDisplay.RescaleTransferFunctionToDataRange()


    renderView1.Update()
    renderView1.ResetCamera()
    renderView1.ViewSize = geometry
    renderView1.CameraPosition = [0.0005945160428284565, 1.959300980464672e-13, -1.3552527156068805e-20]
    renderView1.CameraFocalPoint = [-2.549999918319142e-05, 1.959300980464672e-13, -1.3552527156068805e-20]
    renderView1.CameraViewUp = [0.0, 0.0, 1.0]
    renderView1.CameraParallelScale = 0.00016047195994169908
    renderView1.ResetCamera()

    for colorVar in colorVars:
        print("Animating", colorVar )


        ColorBy(projectionDisplay, ('POINTS', colorVar))
        projectionDisplay.RescaleTransferFunctionToDataRange()

        wLUT = GetColorTransferFunction(colorVar)
        wPWF = GetOpacityTransferFunction(colorVar)
        HideScalarBarIfNotNeeded(wLUT, renderView1)

        ## NOTE: For color presets.
        wLUT.ApplyPreset('Cool to Warm (Extended)', True)

        renderView1.Update()
        UpdateScalarBars()

        # renderView1.ViewSize = [1750, 1300]
        # SaveAnimation(colorVar + '.png', renderView1, ImageResolution=[1750, 1300], TransparentBackground=1, SuffixFormat='.%04d')

        SaveAnimation(colorVar + '.png', renderView1, ImageResolution=geometry, TransparentBackground=1, SuffixFormat='.%04d')


def main():

    ap = argparse.ArgumentParser()

    ap.add_argument("-p", "--projectionType", required=False, default='clip', help="projection type: clip | slice")
    ap.add_argument("-s", "--snapshot", required=False, action='store_true', help="run snapshotter")
    ap.add_argument("-d", "--distribute", required=False, help="Apply d3 filter and save data")
    ap.add_argument("-m", "--mass-defect", required=False, action='store_true', help="Integrate and find mass_defect curve")

    ap.add_argument("-c", "--colorVars", required=False, nargs='*', help="color map variable")
    ap.add_argument("-g", "--geometry", required=False, nargs=2, type=int, default=[1750, 1300], help="Animation geometry size")
    ap.add_argument("-nsb", "--no-scalar-bar", required=False, action='store_true', default=False, help="Disable scalar bar visibility")
    ap.add_argument("-nca", "--no-coordinate-axis", required=False, action='store_true', default=False, help="Disable coordinate axis visibility")

    ap.add_argument("-f", "--filetype", required=False, default='pvtu', choices=['xdmf', 'vtu', 'vtk', 'pvtu'], help="filetype: xdmf | vtu | vtk | pvtu")
    ap.add_argument("-w", "--writer", required=False, choices=['xdmf', 'vtu', 'vtk', 'pvtu'], help="writer: xdmf | vtu | vtk | pvtu")

    ap.add_argument("FILES", nargs='*', help="files..")

    args = vars(ap.parse_args())

    # print(args['geometry'])
    # for x in args['geometry']:
    #     print(type(x))

    if len(args['FILES']) == 0:
        filetype = args['filetype']

        args['FILES'] = sorted([ file for file in os.listdir(".") if file.endswith(filetype) ], key=lambda x: int(x.split('.')[0].split('_')[-1]))

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
        reader = LegacyVTKReader(FileName=args['FILES'])
    else:
        print("Unsupported File Format!")
        sys.exit(-1)

    if args['distribute']:
        d3 = D3(Input=reader)
        SaveData(args['distribute'] + '.pvtu', proxy=d3,
                Writealltimestepsasfileseries=1)
        sys.exit(0)

    if args['snapshot']:
        snapshot(reader,
                projectionType = args['projectionType'],
                colorVars = args['colorVars'],
                scalarBarVisible = not args['no_scalar_bar'],
                geometry = args['geometry'],
                axisVisible = not args['no_coordinate_axis']
                )

    if args['mass_defect']:
        mass_defect(reader);

    if args['writer']:
        writer=None

        if args['writer'] == 'xdmf':
            writer = XDMFWriter(FileNames=args['FILES'])
        elif args['writer'] == 'vtu':
            writer = XMLUnstructuredGridWriter(FileName=args['FILES'])
            writer.Writetimestepsasfileseries=1
        elif args['writer'] == 'pvtu':
            writer = XMLPUnstructuredGridWriter(FileName=args['FILES'])
            writer.Writealltimestepsasfileseries=1
        elif args['writer'] == 'vtk':
            writer = LegacyVTKWriter(FileName=args['FILES'])
        else:
            print("Reader Unspecified!")
            sys.exit(-1)

        writer.Timestepstride = 1
        writer.FileName = 'script-output'
        writer.UpdatePipeline()


if __name__ == "__main__":
    main()
