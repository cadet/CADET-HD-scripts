#!/bin/env -S pvpython --force-offscreen-rendering

## https://www.paraview.org/Wiki/ParaView_and_Python#Control_the_camera

"""
A Paraview script to aid the postprocessing of data. I wrote it mainly to generate images out of output files.

Usage by examples:
    > ./vis.py -h
    > pvpython --force-offscreen-rendering /path/to/vis.py **/*pvtu -a -c scalar_2 scalar_3 -p slice  ## for listed pvtu files, output scalar_2 and scalar_3 slice data
    > mpirun -np 4 pvbatch /path/to/vis.py -f pvtu -a                   ## looks for pvtu files in current folder, clips them, animates all scalars,
    > mpirun -np 4 pvbatch /path/to/vis.py **/*pvtu -d out-distributed  ## for listed pvtu files, apply d3 filter and write to out-distributed
    > mpirun -np 4 pvbatch /path/to/vis.py example-case_{150..300}.pvtu -a  ## Uses bash's sequence expansion

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
## TODO: fix timestep data to be handled in from pvtu file? (in mixd2pvtu)
## TODO: rotated views [Good for image resolution in long columns]
## DONE: After updatescalarbars() they are always visible even with -nsb
## TODO: Allow modularity/composability of functions in operating modes

import argparse
from matplotlib import pyplot as plt
from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()
from vtkmodules.numpy_interface import dataset_adapter as dsa
import vtk.util.numpy_support as ns
from math import asin,sqrt,pi
import numpy as np
import pickle
import struct
import sys
import csv
import os

# import h5py ##NOTE: Is perhaps better as a standard format, but issue with the version of paraview compiled locally at time of scripting. So working with bin files.

def MIXDReader():
    pass

def appendToBin(arr, filename, f):
    with(open(filename, 'ab')) as output:
        for i in arr:
            output.write(struct.pack(f, i))

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

def csvWriter(filename, x, y):
    with open(filename, 'w') as f:
        writer = csv.writer(f)
        writer.writerows(zip(x, y))

## TODO: Test written bin in mpi mode
def bead_loading(reader, **kwargs):
    colorVars = kwargs.get('colorVars', reader.PointArrayStatus) or reader.PointArrayStatus
    files = kwargs.get('files')

    try:
        timeKeeper = GetTimeKeeper()
        nts = len(reader.TimestepValues)
    except:
        nts = 1
        pass

    if nts == 0:
        nts = 1
    ncv = len(colorVars)

    print("nts:", nts)
    print("ncv:", ncv)

    renderView1 = GetActiveViewOrCreate('RenderView')
    # display = Show(reader, renderView1)
    # Hide(reader, renderView1)

    connectivity = Connectivity(Input=reader)
    connectivityDisplay = Show(connectivity, renderView1)
    Hide(connectivity, renderView1)

    # NOTE: Threshold  range will be (0, n) where n is number of beads.
    # Typically, the interstitial domain is the last, n+1th region.
    # Here, we ignore the interstitial region by setting nbeads = n, and not n+1.
    nbeads = int(connectivity.PointData.GetArray("RegionId").GetRange()[1])
    print("Number of Objects:", nbeads)

    appendToBin([nts, nbeads, ncv],'bead_loading.inf', '=i')
    dataArr = np.zeros((nts, nbeads, ncv))
    coordArr = np.zeros((nbeads,4))


    for timestep in range(nts):
        timeKeeper.Time = timestep
        # print("Processing timestep: ", timestep, end="\r")

        for index in range(nbeads):

            print("Processing timestep: {timestep:3d} | bead: {index:5d}".format(timestep=timestep, index=index), end="\r")
            threshold = Threshold(Input=connectivity)
            threshold.ThresholdRange = [index, index]
            thresholdDisplay = Show(threshold, renderView1)
            # threshold.UpdatePipeline()

            if timestep == 0:
                (xmin,xmax,ymin,ymax,zmin,zmax) = GetActiveSource().GetDataInformation().GetBounds()
                # print("Threshold coordinate bounds:",xmin,xmax,ymin,ymax,zmin,zmax)
                x = (xmax + xmin)/2
                y = (ymax + ymin)/2
                z = (zmax + zmin)/2
                r = (xmax - xmin + ymax - ymin + zmax - zmin)/6
                # print("xyzr:",x, y, z, r)
                # coordArr[index,:] = np.array([x, y, z, r])
                appendToBin([x,y,z,r],'bead_loading.xyzr', '=d')

            integrated = IntegrateVariables(Input=threshold)
            intdata = servermanager.Fetch(integrated)
            intdata = dsa.WrapDataObject(intdata)

            values = []
            for colorVar in colorVars:
                value = intdata.PointData[colorVar]
                value = ns.vtk_to_numpy(value)
                values.append(value[0])

            dataArr[timestep,index,:] = np.array(values)
            Hide(threshold, renderView1)

            Delete(integrated)
            Delete(thresholdDisplay)
            Delete(threshold)

        # TODO: this only works with one scalar currently, which is okay for now
        appendToBin(dataArr[timestep,:,:], 'ts_' + str(timestep) + '.dat', "=d")


def calc_beads_loading(connectivity, colorVars, dataArr, nbeads):
    for index in range(nbeads):

        threshold = Threshold(Input=connectivity)
        threshold.ThresholdRange = [index, index]
        thresholdDisplay = Show(threshold, renderView1)
        # threshold.UpdatePipeline()

        integrated = IntegrateVariables(Input=threshold)
        intdata = servermanager.Fetch(integrated)
        intdata = dsa.WrapDataObject(intdata)

        values = []
        for colorVar in colorVars:
            value = intdata.PointData[colorVar]
            value = ns.vtk_to_numpy(value)
            values.append(value[0])

        # dataArr[timestep,index,:] = np.array(values)
        appendToBin(values, 'bead_' + str(index) + '.dat', "=d")

        Delete(threshold)
        del(threshold)

def mass_flux(reader, **kwargs):
    scalarBarVisible = kwargs.get('scalarBarVisible', True)
    geometry = kwargs.get('geometry', [1750, 1300])
    axisVisible = kwargs.get('axisVisible', True)
    colorVars = kwargs.get('colorVars', reader.PointArrayStatus) or reader.PointArrayStatus
    nSlice = kwargs.get('nSlice', 1) or 1

    renderView1 = GetActiveViewOrCreate('RenderView')
    display = Show(reader, renderView1)
    # display.Representation = 'Surface With Edges'

    (xmin,xmax,ymin,ymax,zmin,zmax) = GetActiveSource().GetDataInformation().GetBounds()
    print(xmin,xmax,ymin,ymax,zmin,zmax)

    Hide(reader, renderView1)

    ## NOTE: Only takes takes one color
    colorVar = colorVars[0]
    flowrate = []
    zs = []

    count = 0

    try:
        ## Set timestep to last timestep (last file in series)
        timeKeeper = GetTimeKeeper()
        timeKeeper.Time = reader.TimestepValues[-1]
    except:
        pass

    for zpos in np.linspace(zmin, zmax, nSlice):

        count = count + 1
        print("Loop: ", count, zpos)
        projection = Slice(Input=reader)
        Hide3DWidgets(proxy=projection.SliceType)

        projection.SliceType = 'Plane'
        projection.HyperTreeGridSlicer = 'Plane'
        projection.SliceOffsetValues = [0.0]

        projection.SliceType.Origin = [0.0, 0.0, zpos]
        projection.SliceType.Normal = [0.0, 0.0, -1.0]
        projection.UpdatePipeline()

        ## {{{
        # projectionDisplay = Show(projection, renderView1)
        # projectionDisplay.Representation = 'Surface'
        # # projectionDisplay.Representation = 'Surface With Edges'
        # renderView1.OrientationAxesVisibility = int(axisVisible)
        # projectionDisplay.RescaleTransferFunctionToDataRange()

        # renderView1.Update()
        # renderView1.ViewSize = geometry
        # renderView1.ResetCamera()

        # ColorBy(projectionDisplay, ('POINTS', colorVar))
        # projectionDisplay.RescaleTransferFunctionToDataRange()

        # wLUT = GetColorTransferFunction(colorVar)
        # wPWF = GetOpacityTransferFunction(colorVar)
        # HideScalarBarIfNotNeeded(wLUT, renderView1)

        # ## NOTE: For color presets.
        # wLUT.ApplyPreset('Cool to Warm (Extended)', True)

        # renderView1.Update()
        # UpdateScalarBars()

        # projectionDisplay.SetScalarBarVisibility(renderView1, scalarBarVisible)

        # SaveScreenshot(colorVar + str(count)+'.png', renderView1, ImageResolution=[1750, 1300], TransparentBackground=0)

        # Hide(projection, renderView1)
        ## }}}

        integrated = IntegrateVariables(Input=projection)
        intdata = servermanager.Fetch(integrated)
        intdata = dsa.WrapDataObject(intdata)
        value = intdata.PointData[colorVar]
        try:
            value = ns.vtk_to_numpy(value)
            flowrate.append(value[0])
            zs.append(zpos)
        except:
            pass

    print(zs)
    print(flowrate)
    csvWriter('massFlux.csv', zs, flowrate)
    plt.figure()
    plt.plot(zs, flowrate)
    plt.savefig('plot.pdf')
    # plt.show()

def snapshot(reader, **kwargs):
    scalarBarVisible = kwargs.get('scalarBarVisible', True)
    axisVisible = kwargs.get('axisVisible', True)
    geometry = kwargs.get('geometry')
    zoom = kwargs.get('zoom')

    renderView1 = GetActiveViewOrCreate('RenderView')
    display = Show(reader, renderView1)
    display.Representation = 'Surface'
    # display.Representation = 'Surface With Edges'
    renderView1.OrientationAxesVisibility = int(axisVisible)
    # display.RescaleTransferFunctionToDataRange()

    renderView1.Update()
    renderView1.ResetCamera()
    renderView1.ViewSize = geometry

    setCameraOrientation()

    # renderView1.CameraPosition = [0.0005945160428284565, 1.959300980464672e-13, -1.3552527156068805e-20]
    # renderView1.CameraFocalPoint = [-2.549999918319142e-05, 1.959300980464672e-13, -1.3552527156068805e-20]
    # renderView1.CameraViewUp = [0.0, 0.0, 1.0]
    # renderView1.CameraParallelScale = 0.00016047195994169908
    # renderView1.ResetCamera()

    cam = GetActiveCamera()
    cx,cy,cz = cam.GetPosition()
    cam.SetPosition(cx/zoom[0], cy/zoom[0], cz/zoom[0])
    Render()

    timeKeeper1 = GetTimeKeeper()
    for ts in reader.TimestepValues:
        timeKeeper1.Time = ts
        SaveScreenshot(str(ts)+'.png', renderView1, ImageResolution=[1750, 1300], TransparentBackground=1)

def animate(reader, **kwargs):
    projectionType = kwargs.get('projectionType', 'clip')
    scalarBarVisible = kwargs.get('scalarBarVisible', True)
    axisVisible = kwargs.get('axisVisible', True)
    geometry = kwargs.get('geometry')
    colorVars = kwargs.get('colorVars', reader.PointArrayStatus) or reader.PointArrayStatus

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
    renderView1.OrientationAxesVisibility = int(axisVisible)
    projectionDisplay.RescaleTransferFunctionToDataRange()
    renderView1.ViewSize = geometry
    renderView1.Update()

    setCameraOrientation()

    for colorVar in colorVars:
        print("Animating", colorVar )

        ColorBy(projectionDisplay, ('POINTS', colorVar))
        projectionDisplay.RescaleTransferFunctionToDataRange()

        wLUT = GetColorTransferFunction(colorVar)
        wPWF = GetOpacityTransferFunction(colorVar)
        HideScalarBarIfNotNeeded(wLUT, renderView1)

        ## NOTE: For color presets.
        wLUT.ApplyPreset('Rainbow Uniform', True)

        renderView1.Update()
        UpdateScalarBars()

        projectionDisplay.SetScalarBarVisibility(renderView1, scalarBarVisible)
        SaveAnimation(colorVar + '.png', renderView1, ImageResolution=geometry, TransparentBackground=1, SuffixFormat='.%04d')


def radial_shell_integrate(reader, **kwargs):

    # projectionType = kwargs.get('projectionType', 'clip')
    scalarBarVisible = kwargs.get('scalarBarVisible', True)
    axisVisible = kwargs.get('axisVisible', True)
    geometry = kwargs.get('geometry')
    colorVars = kwargs.get('colorVars', reader.PointArrayStatus) or reader.PointArrayStatus
    zoom = kwargs.get('zoom')
    nRegions = kwargs.get('nRegions')

    # rShells = [0, 1e-5, 2e-5, 3e-5, 4e-5, 5e-5]

    ## Calc bounding box. Requires show
    renderView1 = GetActiveViewOrCreate('RenderView')
    display = Show(reader, renderView1)
    (xmin,xmax,ymin,ymax,zmin,zmax) = GetActiveSource().GetDataInformation().GetBounds()
    Hide(reader, renderView1)

    # nRegions = 3
    nShells = nRegions + 1 #Including r = 0
    rShells = []

    R = (xmax - xmin + ymax - ymin)/4
    print("R:", R)

    shellType = 'EQUIDISTANT'
    # shellType = 'EQUIVOLUME'
    if shellType == 'EQUIVOLUME':
        for n in range(nShells):
            rShells.append(R * sqrt(n/nRegions))
    elif shellType == 'EQUIDISTANT':
        for n in range(nShells):
            rShells.append(R * (n/nRegions))

    print("rShells:", rShells)

    appended = []
    radAvg = []

    for radIn, radOut in zip(rShells[:-1], rShells[1:]+rShells[:0]):

        radAvg.append( (radIn + radOut) / 2 )

        clipOuter = Clip(Input=reader)
        clipOuter.ClipType = 'Cylinder'
        clipOuter.ClipType.Axis = [0.0, 0.0, 1.0]
        clipOuter.ClipType.Radius = radOut
        Hide3DWidgets(proxy=clipOuter.ClipType)

        # renderView1 = GetActiveViewOrCreate('RenderView')
        # projectionDisplay = Show(clipOuter, renderView1)
        # projectionDisplay.Representation = 'Surface'
        # # projectionDisplay.Representation = 'Surface With Edges'
        # renderView1.OrientationAxesVisibility = int(axisVisible)
        # projectionDisplay.RescaleTransferFunctionToDataRange()

        clipInner = Clip(Input=clipOuter)
        clipInner.ClipType = 'Cylinder'
        clipInner.ClipType.Axis = [0.0, 0.0, 1.0]
        clipInner.ClipType.Radius = radIn
        clipInner.Invert = 0

        cellSize1 = CellSize(Input=clipInner)
        cellSize1.ComputeVolume = 1
        cellSize1.ComputeSum = 1

        volume = servermanager.Fetch(cellSize1)
        volume = dsa.WrapDataObject(volume)
        volume = volume.FieldData['Volume'][0]
        print("VOLUME:", volume)

        integrated = IntegrateVariables(Input=clipInner)
        intdata = servermanager.Fetch(integrated)
        intdata = dsa.WrapDataObject(intdata)

        values = []
        for colorVar in colorVars:
            value = intdata.PointData[colorVar]
            value = ns.vtk_to_numpy(value)
            values.append(value[0]/volume) ## Average of velocity, instead of integ(v.dV)

        print(values)
        appended.extend(values)

    print("Average scalar by radius:", appended)
    csvWriter('radial_shell_integrate', radAvg, appended)


# def chromatogram(reader, **kwargs):
#     extractSurface1 = ExtractSurface(Input=testcase_0pvtu)
#     generateSurfaceNormals1 = GenerateSurfaceNormals(Input=extractSurface1)
#     threshold1 = Threshold(Input=generateSurfaceNormals1)
#     threshold1.ThresholdRange = [1.0, 1.0]
#     threshold1.Scalars = ['POINTS', 'Normals_Z']


def main():

    ap = argparse.ArgumentParser()

    ap.add_argument("-p", "--projectionType", required=False, default='clip', help="projection type: clip | slice")
    ap.add_argument("-a", "--animate", required=False, action='store_true', help="run animator")
    ap.add_argument("-d", "--distribute", required=False, help="Apply d3 filter and save data")
    ap.add_argument("-m", "--mass-flux", required=False, help="Find mass flux at n different slices along the z direction")
    ap.add_argument("-s", "--snapshot", required=False, action='store_true', help="run snapshotter")
    ap.add_argument("-b", "--bead-loading", required=False, action='store_true', help="Output bead loading data")
    ap.add_argument("-r", "--radial-integrate", required=False, help="Cylindrical shell integrate variables")

    ap.add_argument("-c", "--colorVars", required=False, nargs='*', help="color map variable")
    ap.add_argument("-g", "--geometry", required=False, nargs=2, type=int, default=[1750, 1300], help="Animation geometry size")
    ap.add_argument("-z", "--zoom", required=False, nargs=3, type=float, default=[1, 1, 1], help="Zoom factors for snapshot")
    ap.add_argument("-nsb", "--no-scalar-bar", required=False, action='store_true', default=False, help="Disable scalar bar visibility")
    ap.add_argument("-nca", "--no-coordinate-axis", required=False, action='store_true', default=False, help="Disable coordinate axis visibility")
    ap.add_argument("-css", "--cross-section-snapshots", type=int, required=False, help="Run snapshotter for n cross section slices")

    ap.add_argument("-f", "--filetype", required=False, default='pvtu', choices=['xdmf', 'vtu', 'vtk', 'pvtu'], help="filetype: xdmf | vtu | vtk | pvtu")
    ap.add_argument("-w", "--writer", required=False, choices=['xdmf', 'vtu', 'vtk', 'pvtu'], help="writer: xdmf | vtu | vtk | pvtu")

    ap.add_argument("FILES", nargs='*', help="files..")

    args = vars(ap.parse_args())

    # print(args['geometry'])
    # for x in args['geometry']:
    #     print(type(x))

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
        sys.exit(-1)

    if args['distribute']:
        d3 = D3(Input=reader)
        SaveData(args['distribute'] + '.pvtu', proxy=d3,
                Writealltimestepsasfileseries=1)
        sys.exit(0)

    if args['animate']:
        animate(reader,
                projectionType = args['projectionType'],
                colorVars = args['colorVars'],
                scalarBarVisible = not args['no_scalar_bar'],
                geometry = args['geometry'],
                axisVisible = not args['no_coordinate_axis']
                )

    if args['snapshot']:
        snapshot(reader,
                scalarBarVisible = not args['no_scalar_bar'],
                geometry = args['geometry'],
                axisVisible = not args['no_coordinate_axis'],
                zoom = args['zoom']
                )

    if args['mass_flux']:
        mass_flux(reader,
                scalarBarVisible = not args['no_scalar_bar'],
                geometry = args['geometry'],
                axisVisible = not args['no_coordinate_axis'],
                colorVars = args['colorVars'],
                nSlice = args['mass_flux']
                )

    if args['bead_loading']:
        bead_loading(reader,
                colorVars = args['colorVars'],
                files = args['FILES']
                )

    if args['radial_integrate']:
        radial_shell_integrate( reader,
                scalarBarVisible = not args['no_scalar_bar'],
                geometry = args['geometry'],
                axisVisible = not args['no_coordinate_axis'],
                zoom = args['zoom'],
                colorVars = args['colorVars'],
                nRegions = int(args['radial_integrate'])
                )

    if args['cross_section_snapshots']:
        cross_section_snapshots(reader,
                scalarBarVisible = not args['no_scalar_bar'],
                geometry = args['geometry'],
                axisVisible = not args['no_coordinate_axis'],
                zoom = args['zoom'],
                colorVars = args['colorVars'],
                nSlice = args['cross_section_snapshots'],
                files = args['FILES']
                )


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

def cross_section_snapshots(reader, **kwargs):
    scalarBarVisible = kwargs.get('scalarBarVisible', True)
    geometry = kwargs.get('geometry', [1750, 1300])
    axisVisible = kwargs.get('axisVisible', True)
    colorVars = kwargs.get('colorVars', reader.PointArrayStatus) or reader.PointArrayStatus
    nSlice = kwargs.get('nSlice', 1) or 1
    files = kwargs.get('files')


    renderView1 = GetActiveViewOrCreate('RenderView')
    display = Show(reader, renderView1)
    # display.Representation = 'Surface With Edges'

    (xmin,xmax,ymin,ymax,zmin,zmax) = GetActiveSource().GetDataInformation().GetBounds()
    print(xmin,xmax,ymin,ymax,zmin,zmax)

    Hide(reader, renderView1)

    ## NOTE: Only takes takes one color
    colorVar = colorVars[0]
    flowrate = []
    zs = []

    count = 0

    try:
        ## Set timestep to last timestep (last file in series)
        timeKeeper = GetTimeKeeper()
        timeKeeper.Time = reader.TimestepValues[-1]
    except:
        pass

    for zpos in np.linspace(zmin, zmax, nSlice):

        count = count + 1
        print("Loop: ", count, zpos)
        projection = Slice(Input=reader)
        Hide3DWidgets(proxy=projection.SliceType)

        projection.SliceType = 'Plane'
        projection.HyperTreeGridSlicer = 'Plane'
        projection.SliceOffsetValues = [0.0]

        projection.SliceType.Origin = [0.0, 0.0, zpos]
        projection.SliceType.Normal = [0.0, 0.0, -1.0]
        projection.UpdatePipeline()

        projectionDisplay = Show(projection, renderView1)
        projectionDisplay.Representation = 'Surface'
        # projectionDisplay.Representation = 'Surface With Edges'
        renderView1.OrientationAxesVisibility = int(axisVisible)
        projectionDisplay.RescaleTransferFunctionToDataRange()

        renderView1.Update()
        renderView1.ViewSize = geometry
        renderView1.ResetCamera()

        ColorBy(projectionDisplay, ('POINTS', colorVar))
        projectionDisplay.RescaleTransferFunctionToDataRange()

        wLUT = GetColorTransferFunction(colorVar)
        wPWF = GetOpacityTransferFunction(colorVar)
        HideScalarBarIfNotNeeded(wLUT, renderView1)

        ## NOTE: For color presets.
        wLUT.ApplyPreset('Rainbow Uniform', True)

        renderView1.Update()
        UpdateScalarBars()

        projectionDisplay.SetScalarBarVisibility(renderView1, scalarBarVisible)

        SaveScreenshot(colorVar + '_' + str(count)+'.png', renderView1, ImageResolution=[1750, 1300], TransparentBackground=1)

        Hide(projection, renderView1)


def setCameraOrientation():
    camera = GetActiveCamera()
    camera.SetFocalPoint(-1,0,0)
    camera.SetPosition(1,0,0)
    camera.SetViewUp(0,-1,0)
    ResetCamera()

if __name__ == "__main__":
    main()
