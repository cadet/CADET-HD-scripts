#!/bin/env -S pvpython --force-offscreen-rendering

## https://www.paraview.org/Wiki/ParaView_and_Python#Control_the_camera

""" vis.py
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
    > -f filetype is optional. Default is pvtu. This way, running vis.py in a folder lets it detect & sort all files of the given filetype in the current folder unless files are specifically input.
    > The program ALWAYS starts animation outputs with suffix 0000, regardless of the timestep data you've supplied. So be wary of using it twice in the same dir.
    > Assumes that timestep suffix is "_%d" to autodetect files in current directory

"""

## TODO: config.json file
## TODO: Allow modularity/composability of functions in operating modes
## TODO: Automatic filetype detection
## TODO: Easy way to select "Solid Color"


## TODO: get min/max over all timesteps for scalarBar for animations etc
# You could run Temporal Statistics to get Max and Min over entire time range, then run calculator to get Range as (Max-Min).
# Then you could run your filter on this and use the variable from the calculator.

import argparse
from matplotlib import pyplot as plt
# from paraview.simple import * #type:ignore
from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()
from vtkmodules.numpy_interface import dataset_adapter as dsa
import vtk.util.numpy_support as ns #type:ignore
from math import sqrt
import numpy as np
import pickle
import struct
import sys
import csv
import os

# import h5py ##NOTE: Is perhaps better as a standard format, but issue with the version of paraview compiled locally at time of scripting. So working with bin files.


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
def bead_loading(reader, args):
    colorVars = args['colorVars'] or reader.PointArrayStatus
    files = args['FILES']

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

    view = GetActiveViewOrCreate('RenderView')

    connectivity = Connectivity(Input=reader)
    connectivityDisplay = Show(connectivity, view)
    Hide(connectivity, view)

    # NOTE: Threshold  range will be (0, n) where n is number of beads.
    # Typically, the interstitial domain is the last, n+1th region.
    # Here, we ignore the interstitial region by setting nbeads = n, and not n+1.
    nbeads = int(connectivity.PointData.GetArray("RegionId").GetRange()[1])
    print("Number of Objects:", nbeads)

    appendToBin([nts, nbeads, ncv],'bead_loading.inf', '=i')
    dataArr = np.zeros((nts, nbeads, ncv))
    # coordArr = np.zeros((nbeads,4))


    for timestep in range(nts):
        timeKeeper.Time = timestep
        # print("Processing timestep: ", timestep, end="\r")

        for index in range(nbeads):

            print("Processing timestep: {timestep:3d} | bead: {index:5d} | file: {file}".format(timestep=timestep, index=index, file=files[timestep]), end="\r")
            threshold = Threshold(Input=connectivity)
            threshold.ThresholdRange = [index, index]
            thresholdDisplay = Show(threshold, view)
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
            Hide(threshold, view)

            Delete(integrated)
            Delete(thresholdDisplay)
            Delete(threshold)

        # TODO: this only works with one scalar currently, which is okay for now
        # appendToBin(dataArr[timestep,:,:], 'ts_' + str(timestep) + '.dat', "=d")
        appendToBin(dataArr[timestep,:,:], files[timestep].replace('.pvtu', '.dat'), "=d")


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

def mass_flux(reader, args):

    colorVars = args['colorVars'] or reader.PointArrayStatus
    nSlice = args['mass_flux'] or 1

    view = GetActiveViewOrCreate('RenderView')
    display = Show(reader, view)
    display.Representation = args['display_representation']
    # display.Representation = 'Surface With Edges'

    (xmin,xmax,ymin,ymax,zmin,zmax) = GetActiveSource().GetDataInformation().GetBounds()
    print(xmin,xmax,ymin,ymax,zmin,zmax)

    Hide(reader, view)

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

def geometry_snapshot(reader, args):

    geometry    = args['geometry']
    axisVisible = not args['no_coordinate_axis']
    zoom        = args['zoom']
    files       = args['FILES']
    filetype    = args['filetype']

    view = GetActiveViewOrCreate('RenderView')
    connectivity = Connectivity(Input=reader)
    connectivityDisplay = Show(connectivity, view)

    connectivityDisplay.Representation = args['display_representation']
    view.OrientationAxesVisibility = int(axisVisible)


    Hide(connectivity, view)

    # NOTE: Threshold  range will be (0, n) where n is number of beads.
    # Typically, the interstitial domain is the last, n+1th region.
    # Here, we ignore the interstitial region by setting nbeads = n, and not n+1.
    nbeads = int(connectivity.PointData.GetArray("RegionId").GetRange()[1])
    print("Number of Objects:", nbeads)

    ## TODO: just use 0..n as threshold range instead of this nonsense
    # for index in range(nbeads):
    # print("Processing bead: {index}".format(index=index))
    print("Processing beads: {nbeads}".format(nbeads=nbeads))
    threshold = Threshold(Input=connectivity)
    threshold.ThresholdRange = [0, nbeads-1]
    thresholdDisplay = Show(threshold, view)
    ColorBy(thresholdDisplay, None)
    # threshold.UpdatePipeline()
    thresholdDisplay.AmbientColor = [2/255, 61/255, 107/255]
    thresholdDisplay.DiffuseColor = [2/255, 61/255, 107/255]

    print("Processing Column.")
    threshold = Threshold(Input=connectivity)
    threshold.ThresholdRange = [nbeads, nbeads]
    outerShell = Projection(threshold, 'clip')
    outerShellDisplay = Show(outerShell, view)
    ColorBy(outerShellDisplay, None)
    # outerShell.UpdatePipeline()
    outerShellDisplay.Opacity = 0.5
    view.InteractionMode = '2D'

    view.Update()
    view.ResetCamera()
    view.ViewSize = geometry
    setCameraOrientation(zoom) ## NOTE: Zoom doesn't work with InteractionMode = '2D'

    ## NOTE: Only works on the first file provided
    SaveScreenshot(files[0].replace(filetype, 'png'), view, ImageResolution=geometry, TransparentBackground=1)

def snapshot(reader, args):
    colorVars = args['colorVars'] or reader.PointArrayStatus
    scalarBarVisible = not args['no_scalar_bar']
    geometry = args['geometry']
    axisVisible = not args['no_coordinate_axis']
    zoom = args['zoom']
    files = args['FILES']
    filetype = args['filetype']
    projectionType = args['projectionType']


    ## TODO: implement cutting plane coordinates
    projection = Projection(reader, projectionType)
    # Render()

    view = GetActiveViewOrCreate('RenderView')
    # renderView1 = GetActiveView()
    display = Show(projection, view)
    display.Representation = args['display_representation']
    view.OrientationAxesVisibility = int(axisVisible)

    view.Update()
    view.ResetCamera()
    view.ViewSize = geometry

    display.UpdatePipeline()

    setCameraOrientation(zoom)

    timeKeeper1 = GetTimeKeeper()

    for ifile in files:
        try:
            timeKeeper1.Time = reader.TimestepValues[files.index(ifile)]
        except:
            timeKeeper1.Time = 0

        for colorVar in colorVars:
            print("Snapping", colorVar )

            if colorVar == 'None':
                ColorBy(display, None)
            else:
                ColorBy(display, ('POINTS', colorVar))

            display.RescaleTransferFunctionToDataRange()

            wLUT = GetColorTransferFunction(colorVar)
            wPWF = GetOpacityTransferFunction(colorVar)
            HideScalarBarIfNotNeeded(wLUT, view)

            wLUT.ApplyPreset('Rainbow Uniform', True)

            view.Update()
            UpdateScalarBars()

            display.SetScalarBarVisibility(view, scalarBarVisible)
            SaveScreenshot(colorVar + '_' + ifile.replace(filetype, 'png'), view, ImageResolution=geometry, TransparentBackground=1)

def animate(reader, args):
    projectionType = args['projectionType']
    colorVars = args['colorVars'] or reader.PointArrayStatus
    scalarBarVisible = not args['no_scalar_bar']
    geometry = args['geometry']
    axisVisible = not args['no_coordinate_axis']
    zoom = args['zoom']

    animationScene1 = GetAnimationScene()
    timeKeeper1 = GetTimeKeeper()
    animationScene1.UpdateAnimationUsingDataTimeSteps()

    ## TODO: Animate using constant scalarbar range
    ## TODO: Fix animation for one timestep

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

    projection = Projection(reader, projectionType)

    view = GetActiveViewOrCreate('RenderView')
    projectionDisplay = Show(projection, view)
    projectionDisplay.Representation = args['display_representation']
    view.OrientationAxesVisibility = int(axisVisible)
    projectionDisplay.RescaleTransferFunctionToDataRange()
    view.ViewSize = geometry
    view.Update()

    setCameraOrientation(zoom)

    for colorVar in colorVars:
        print("Animating", colorVar )

        if colorVar == 'None':
            ColorBy(projectionDisplay, None)
        else:
            ColorBy(projectionDisplay, ('POINTS', colorVar))

        ## NOTE: Removing this should HELP fix the varying scalar bar range for every frame
        projectionDisplay.RescaleTransferFunctionToDataRange()

        wLUT = GetColorTransferFunction(colorVar)
        wPWF = GetOpacityTransferFunction(colorVar)
        HideScalarBarIfNotNeeded(wLUT, view)

        wLUT.ApplyPreset('Rainbow Uniform', True)

        view.Update()
        UpdateScalarBars()

        projectionDisplay.SetScalarBarVisibility(view, scalarBarVisible)
        SaveAnimation(colorVar + '.png', view, ImageResolution=geometry, TransparentBackground=1, SuffixFormat='.%04d')


def radial_shell_integrate(reader, args):

    colorVars = args['colorVars'] or reader.PointArrayStatus
    nRegions = int(args['radial_integrate'])

    ## Calc bounding box. Requires show
    view = GetActiveViewOrCreate('RenderView')
    display = Show(reader, view)
    (xmin,xmax,ymin,ymax,zmin,zmax) = GetActiveSource().GetDataInformation().GetBounds()
    Hide(reader, view)

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




def cross_section_snapshots(reader, args):
    scalarBarVisible = not args['no_scalar_bar']
    geometry = args['geometry']
    axisVisible = not args['no_coordinate_axis']
    colorVars = args['colorVars'] or reader.PointArrayStatus
    nSlice = args['cross_section_snapshots'] or 1

    view = GetActiveViewOrCreate('RenderView')
    display = Show(reader, view)
    display.Representation = args['display_representation']
    # display.Representation = 'Surface With Edges'

    (xmin,xmax,ymin,ymax,zmin,zmax) = GetActiveSource().GetDataInformation().GetBounds()
    print(xmin,xmax,ymin,ymax,zmin,zmax)

    Hide(reader, view)

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

        projectionDisplay = Show(projection, view)
        projectionDisplay.Representation = args['display_representation']
        # projectionDisplay.Representation = 'Surface With Edges'
        view.OrientationAxesVisibility = int(axisVisible)
        projectionDisplay.RescaleTransferFunctionToDataRange()

        view.Update()
        view.ViewSize = geometry
        view.ResetCamera()

        ColorBy(projectionDisplay, ('POINTS', colorVar))
        projectionDisplay.RescaleTransferFunctionToDataRange()

        wLUT = GetColorTransferFunction(colorVar)
        wPWF = GetOpacityTransferFunction(colorVar)
        HideScalarBarIfNotNeeded(wLUT, view)

        ## NOTE: For color presets.
        wLUT.ApplyPreset('Rainbow Uniform', True)

        view.Update()
        UpdateScalarBars()

        projectionDisplay.SetScalarBarVisibility(view, scalarBarVisible)

        SaveScreenshot(colorVar + '_' + str(count)+'.png', view, ImageResolution=geometry, TransparentBackground=1)

        Hide(projection, view)


def setCameraOrientation(zoom):
    camera = GetActiveCamera()
    camera.SetFocalPoint(-1,0,0)
    camera.SetPosition(1,0,0)
    camera.SetViewUp(0,-1,0)
    ResetCamera()

    camera.Dolly(zoom)
    Render()

    # cx,cy,cz = camera.GetPosition()
    # camera.SetPosition(cx/zoom[0], cy/zoom[0], cz/zoom[0])
    # Render()

## Ideally Projection(inputView, projectionType, direction, depth fraction)
def Projection(inputView, projectionType):

    renderViewForProjection = GetActiveViewOrCreate('RenderView')
    display = Show(inputView, renderViewForProjection)
    (xmin,xmax,ymin,ymax,zmin,zmax) = GetActiveSource().GetDataInformation().GetBounds()
    Hide(inputView, renderViewForProjection)
    center = [ (xmax+xmin)/2, (ymax+ymin)/2, (zmax+zmin)/2,]

    projection = None
    if projectionType == 'clip':
        projection = Clip(Input=inputView)
        projection.ClipType.Origin = center
        Hide3DWidgets(proxy=projection.ClipType)
    elif projectionType == 'slice':
        projection = Slice(Input=inputView)
        projection.SliceType.Origin = center
        Hide3DWidgets(proxy=projection.SliceType)
    else:
        projection = inputView
    return projection

def chromatogram(reader, args):

    colorVars = args['colorVars'] or reader.PointArrayStatus

    timeKeeper = GetTimeKeeper()
    nts = len(reader.TimestepValues)

    renderView1 = GetActiveViewOrCreate('RenderView')
    surfaces = ExtractSurface(Input=reader)
    surfaceNormals = GenerateSurfaceNormals(Input=surfaces)

    threshold = Threshold(Input=surfaceNormals)
    threshold.ThresholdRange = [1.0, 1.0]
    threshold.Scalars = ['POINTS', 'Normals_Z']

    cellSize1 = CellSize(Input=threshold)
    cellSize1.ComputeArea= 1
    cellSize1.ComputeSum = 1

    area = servermanager.Fetch(cellSize1)
    area = dsa.WrapDataObject(area)
    area = area.FieldData['Area'][0]
    print("AREA:", area)
    print("Note: calculating integral for only first scalar: scalar_0")

    values = []

    for timestep in range(nts):

        timeKeeper.Time = timestep
        threshold.UpdatePipeline(reader.TimestepValues[timestep])

        integrated = IntegrateVariables(Input=threshold)
        intdata = servermanager.Fetch(integrated)
        intdata = dsa.WrapDataObject(intdata)

        # for colorVar in colorVars:
        value = intdata.PointData[colorVars[0]]
        value = ns.vtk_to_numpy(value)
        values.append(value[0]/area)  ## Average of c, instead of integ(c.dV)

    print(values)

    csvWriter('chromatogram.csv', reader.TimestepValues, values )

def volumeIntegral(reader, args):
    colorVars = args['colorVars'] or reader.PointArrayStatus

    timeKeeper = GetTimeKeeper()
    nts = len(reader.TimestepValues)

    view = GetActiveViewOrCreate('RenderView')

    # cellSize1 = CellSize(Input=reader)
    # cellSize1.ComputeVolume= 1
    # cellSize1.ComputeSum = 1
    # volume = servermanager.Fetch(cellSize1)
    # volume = dsa.WrapDataObject(volume)
    # volume = volume.FieldData['Volume'][0]
    # print("volume:", volume)

    intOutput = {}
    for colorVar in colorVars:
        intOutput.update({colorVar: []})
        print(type(intOutput[colorVar]))

    for timestep in range(nts):

        timeKeeper.Time = timestep
        # print("Processing timestep: ", timestep, end="\r")
        reader.UpdatePipeline(reader.TimestepValues[timestep])

        integrated = IntegrateVariables(Input=reader)
        intdata = servermanager.Fetch(integrated)
        intdata = dsa.WrapDataObject(intdata)

        volume = intdata.CellData['Volume'][0]

        for colorVar in colorVars:
            value = intdata.PointData[colorVar]
            value = ns.vtk_to_numpy(value)[0]/volume
            intOutput[colorVar].append(value)

        Delete(integrated)

    print(intOutput)
    for colorVar in colorVars:
        csvWriter(colorVar + '.integrated.csv', reader.TimestepValues, intOutput[colorVar])

def main():

    ap = argparse.ArgumentParser()

    ap.add_argument("-s", "--snapshot", required=False, action='store_true', help="run snapshotter")
    ap.add_argument("-a", "--animate", required=False, action='store_true', help="run animator")
    ap.add_argument("-d", "--distribute", required=False, help="Apply d3 filter and save data")
    ap.add_argument("-m", "--mass-flux", required=False, help="Find mass flux at n different slices along the z direction")
    ap.add_argument("-b", "--bead-loading", required=False, action='store_true', help="Output bead loading data")
    ap.add_argument("-r", "--radial-integrate", required=False, help="Cylindrical shell integrate variables")
    ap.add_argument("-css", "--cross-section-snapshots", type=int, required=False, help="Run snapshotter for n cross section slices")
    ap.add_argument("-gs", "--geometry-snapshot", required=False, action='store_true', help="Run snapshotter for n cross section slices")
    ap.add_argument("-p", "--projectionType", required=False, help="projection type: clip | slice")
    ap.add_argument("-cg", "--chromatogram", required=False, action='store_true', help="extract chromatogram from surface with N=(0,0,1)")
    ap.add_argument("-vi", "--volumeIntegral", required=False, action='store_true', help="Calculate volume integral for scalar values in the whole domain")

    ap.add_argument("-c", "--colorVars", required=False, nargs='*', help="color map variable")
    ap.add_argument("-g", "--geometry", required=False, nargs=2, type=int, default=[1750, 1300], help="Animation geometry size")
    ap.add_argument("-z", "--zoom", required=False, type=float, default=1, help="Zoom (camera.dolly) value for view")
    ap.add_argument("-dr", "--display-representation", required=False, default='Surface', choices=['Surface With Edges', 'Surface'], help="Display representation")
    ap.add_argument("-nsb", "--no-scalar-bar", required=False, action='store_true', default=False, help="Disable scalar bar visibility")
    ap.add_argument("-nca", "--no-coordinate-axis", required=False, action='store_true', default=False, help="Disable coordinate axis visibility")

    ap.add_argument("-f", "--filetype", required=False, default='pvtu', choices=['xdmf', 'vtu', 'vtk', 'pvtu'], help="filetype: xdmf | vtu | vtk | pvtu")
    ap.add_argument("-w", "--writer", required=False, choices=['xdmf', 'vtu', 'vtk', 'pvtu'], help="writer: xdmf | vtu | vtk | pvtu")

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
        sys.exit(-1)

    if args['distribute']:
        d3 = D3(Input=reader)
        SaveData(args['distribute'] + '.pvtu', proxy=d3,
                Writealltimestepsasfileseries=1)
        sys.exit(0)

    if args['animate']:
        animate(reader, args)

    if args['snapshot']:
        snapshot(reader, args)

    if args['mass_flux']:
        mass_flux(reader, args)

    if args['bead_loading']:
        bead_loading(reader, args)

    if args['radial_integrate']:
        radial_shell_integrate( reader, args)

    if args['cross_section_snapshots']:
        cross_section_snapshots(reader, args)

    if args['geometry_snapshot']:
        geometry_snapshot(reader, args)

    if args['chromatogram']:
        chromatogram(reader, args)

    if args['volumeIntegral']:
        volumeIntegral(reader, args)

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
