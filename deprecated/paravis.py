#!/bin/env -S pvpython --force-offscreen-rendering

## https://www.paraview.org/Wiki/ParaView_and_Python#Control_the_camera

""" The successor script to vis.py
NOTE: THIS SCRIPT IS DEPRECATED. USE PARAVISION INSTEAD.

    @ideal usage:
        $ paravis.py --pipeline project screenshot --project clip Plane 0.1 z -s None -v z y

"""

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

def csvWriter(filename, x, y):
    with open(filename, 'w') as f:
        writer = csv.writer(f)
        writer.writerows(zip(x, y))

def arr_to_bin(arr, filename, dataformat):
    # NOTE: see https://docs.python.org/3/library/struct.html
    with(open(filename, 'wb')) as output:
        for i in arr:
            output.write(struct.pack(dataformat, i))

def arr_to_bin_unpacked(arr, filename, vartype, mode:str='w'):
    ## vartype = d, f etc
    ## NOTE: uses native by default endianness
    ## Probably faster than the default one
    with(open(filename, mode+'b')) as output:
        output.write(struct.pack(str(len(arr)) + vartype, *arr))

def reference_axisNormal(input:str):
    if input.lower() == "x":
        return [1, 0, 0]
    elif input.lower() == "y":
        return [0, 1, 0]
    elif input.lower() == "z":
        return [0, 0, 1]

# def project(object, projectionType, geometry='Plane', origin=None, normal=[1,0,0]):
def project(object, args):

    projectionType = args['project'][0]
    geometry = args['project'][1]
    origin = args['project'][2]
    normal = args['project'][3]

    origin, normal = default_origin_normal(object, origin, normal)

    # projectionView = GetActiveViewOrCreate('RenderView')
    # display = Show(object, projectionView)
    # (xmin,xmax,ymin,ymax,zmin,zmax) = GetActiveSource().GetDataInformation().GetBounds()
    # Hide(object, projectionView)
    # center = [ (xmax+xmin)/2, (ymax+ymin)/2, (zmax+zmin)/2,]
    # if not origin:
    #     origin=center
    # else:


    projection = None
    if projectionType.lower() == 'clip':
        projection = Clip(Input=object)
        projection.ClipType = geometry
        # projection.HyperTreeGridSlicer = geometry
        projection.ClipType.Origin = origin
        projection.ClipType.Normal = normal
        Hide3DWidgets(proxy=projection.ClipType)
    elif projectionType.lower() == 'slice':
        projection = Slice(Input=object)
        projection.SliceType = geometry
        projection.HyperTreeGridSlicer = geometry
        projection.SliceType.Origin = origin
        projection.SliceType.Normal = normal
        Hide3DWidgets(proxy=projection.SliceType)
    else:
        projection = object

    projection.UpdatePipeline()

    return projection

def GRM2D(object, args):
    ## Split into axial columns
    ## Split into cylindrical columns
    ## Integrate

    view = GetActiveViewOrCreate('RenderView')
    display = Show(object, view)
    display.Representation = args['display_representation']
    (xmin,xmax,ymin,ymax,zmin,zmax) = GetActiveSource().GetDataInformation().GetBounds()
    radius = ((xmax-xmin) + (ymax-ymin)) / 4
    length = zmax - zmin
    print("Length: {}, Radius: {}".format(length, radius))

    # nCol = 10              # Number of axial regions
    # nRad = 5               # Number of radial regions
    nCol = args['grm_2d'][0]
    nRad = args['grm_2d'][1]

    # dx = length/nCol
    # dr = radius/nRad

    colEdges = np.linspace(zmin, zmax, nCol+1)
    radEdges = np.linspace(0,radius,nRad+1) if args['shelltype'] == 'EQUIDISTANT' else list(x/nRad * radius for x in range(nRad+1))


    # Hide(object, view)
    nColEdgeFractions = np.linspace(0,1,nCol+1)
    nRadEdgeFractions = np.linspace(0,1,nRad+1)

    ## TODO: Make these function arguments
    timeKeeper = GetTimeKeeper()
    timeArray = object.TimestepValues
    nts = len(timeArray) or 1

    # ## NOTE: Object must be reader
    # timeArray = object.TimestepValues

    ## Output vector. Should contain (nts X nCol X nRad X nScalar)
    grm2d_output = []

    grm2d_output_filename = 'grm2d_output.bin'

    ## Hack to remove previous file
    arr_to_bin([], grm2d_output_filename, 'd')

    for timestep in range(nts):

        timeKeeper.Time = timestep
        # object.UpdatePipeline(reader.TimestepValues[timestep])
        object.UpdatePipeline(timeArray[timestep])
        # object.UpdatePipeline()

        grm2d_timestep_output = []

        print("--> TS: {}".format(timestep))

        # for leftEdge, rightEdge in zip(colEdges[:-1], colEdges[1:]):
        for leftEdge, rightEdge in zip(nColEdgeFractions[:-1], nColEdgeFractions[1:]):
            print("  |--> Col: {}/{}".format(np.where(nColEdgeFractions == rightEdge)[0][0],nCol))
            SetActiveSource(object)
            print('[{}, {}]'.format(leftEdge, rightEdge))

            # clipLeftArgs = { 'project' : ['clip', 'Plane', leftEdge , '-z'] }
            # clipRightArgs = { 'project' : ['clip', 'Plane', rightEdge, '+z'] }
            # clipLeft = project(object, clipLeftArgs)
            # clipRight = project(clipLeft, clipRightArgs)
            ## screenshot(clipRight, args, suffix=str(leftEdge) + '_')

            clipBox = Clip(Input=object)
            clipBox.ClipType = 'Box'
            clipBox.Exact = 1
            clipBox.ClipType.UseReferenceBounds = 1
            clipBox.ClipType.Bounds = [0.0, 1.0, 0.0, 1.0, leftEdge, rightEdge]
            # clipBox.UpdatePipeline()

            # Hide(object, view)
            # screenshot(clipBox, args, suffix=str(leftEdge) + '_')

            radAvg = []

            for radIn, radOut in zip(radEdges[:-1], radEdges[1:]):
                radAvg.append( (radIn + radOut) / 2 )
                # print('--> [{}, {}]: {}'.format(radIn, radOut, (radIn+radOut)/2))

                print('    |--> Rad: {}/{}'.format(np.where(radEdges == radOut)[0][0], nRad))

                # clipOuter = Clip(Input=clipRight)
                clipOuter = Clip(Input=clipBox)
                clipOuter.ClipType = 'Cylinder'
                clipOuter.ClipType.Axis = [0.0, 0.0, 1.0]
                clipOuter.ClipType.Radius = radOut
                Hide3DWidgets(proxy=clipOuter.ClipType)

                # renderView1 = GetActiveViewOrCreate('RenderView')
                # projectionDisplay = Show(clipOuter, renderView1)
                # projectionDisplay.Representation = 'Surface'
                # # projectionDisplay.Representation = 'Surface With Edges'
                # renderView1.OrientationAxesVisibility = int(args['show_axis'])
                # projectionDisplay.RescaleTransferFunctionToDataRange()

                clipInner = Clip(Input=clipOuter)
                clipInner.ClipType = 'Cylinder'
                clipInner.ClipType.Axis = [0.0, 0.0, 1.0]
                clipInner.ClipType.Radius = radIn
                clipInner.Invert = 0

                # renderView1 = GetActiveViewOrCreate('RenderView')
                # projectionDisplay = Show(clipInner, renderView1)
                # projectionDisplay.Representation = 'Surface'
                # # projectionDisplay.Representation = 'Surface With Edges'
                # renderView1.OrientationAxesVisibility = int(args['show_axis'])
                # projectionDisplay.RescaleTransferFunctionToDataRange()

                # screenshot(clipInner, args, suffix=str(radIn) + "_")

                integrated_scalars = integrate(clipInner, args['scalars'], normalize='Volume')
                # print('---->', integrated_scalars[0])
                grm2d_timestep_output.extend(integrated_scalars[0])

                Delete(clipInner)
                Delete(clipOuter)

            # Delete(clipLeft)
            # Delete(clipRight)
            Delete(clipBox)

        grm2d_output.extend(grm2d_timestep_output)
        arr_to_bin_unpacked(grm2d_timestep_output, 'grm2d_appended.bin', 'd', mode='a')

    # print(grm2d_output)
    ## NOTE: Uncomment one of the below to save from RAM to disk
    # arr_to_bin(grm2d_output, grm2d_output_filename, 'd')
    # arr_to_bin_unpacked(grm2d_output, grm2d_output_filename, 'd')
    print("DONE!")


def screenshot(object, args, suffix=''):
    view = GetActiveViewOrCreate('RenderView')

    display = Show(object, view)
    display.Representation = args['display_representation']
    view.OrientationAxesVisibility = args['show_axis']

    # view.Update()
    # view.ResetCamera()
    # view.ViewSize = args['geometry']

    view_handler(args['view'], args['zoom'])

    for scalar in args['scalars']:
        print("Snapping", scalar )

        display = Show(object, view)
        display.Representation = args['display_representation']
        view.OrientationAxesVisibility = args['show_axis']
        display.RescaleTransferFunctionToDataRange()
        display.UpdatePipeline()

        view.Update()
        view.ViewSize = args['geometry']
        view.ResetCamera()
        UpdateScalarBars()

        if scalar == 'None':
            ColorBy(display, None)
        else:
            ColorBy(display, ('POINTS', scalar))


        wLUT = GetColorTransferFunction(scalar)
        wPWF = GetOpacityTransferFunction(scalar)
        HideScalarBarIfNotNeeded(wLUT, view)

        wLUT.ApplyPreset('Rainbow Uniform', True)

        view.Update()
        UpdateScalarBars()

        display.SetScalarBarVisibility(view, args['show_scalar_bar'])
        SaveScreenshot('screenshot_' + scalar + '_' + suffix + args['FILES'][0].replace(args['filetype'], 'png'), view, ImageResolution=args['geometry'], TransparentBackground=1)
        Hide(display, view)


def integrate(object, vars, normalize=None, timeArray=[]):
    ## normalize= "Volume" or "Area" or None
    choices = ['Volume', 'Area']

    timeKeeper = GetTimeKeeper()
    nts = len(timeArray) or 1

    # GetActiveViewOrCreate('RenderView')

    VolumeOrArea = []

    integrated_over_time = []
    for timestep in range(nts):
        try:
            timeKeeper.Time = timestep
            object.UpdatePipeline(timeArray[timestep])
        except IndexError:
            pass

        print(f"Integrating timestep: {timestep}")

        integrated = IntegrateVariables(Input=object)
        intdata = servermanager.Fetch(integrated)
        intdata = dsa.WrapDataObject(intdata)

        if not intdata:
            raise(AssertionError)

        volume=1
        if normalize in choices:
            if normalize in intdata.CellData.keys():
                volume = intdata.CellData[normalize][0]
            else:
                print("".join(["Cannot normalize by ", normalize, ". No such CellData!"]))

        # print("{key}: {value}".format(key=normalize, value=volume))
        # VolumeOrArea.append(volume)
        # print(VolumeOrArea)

        integrated_scalars = []
        for var in vars:
            value = intdata.PointData[var]
            value = ns.vtk_to_numpy(value)
            integrated_scalars.append(value[0]/volume)  ## Average of c, instead of integ(c.dV)

        integrated_over_time.append(integrated_scalars)

        Delete(integrated)

    return integrated_over_time


def get_cross_sections(reader, nSlice=1):
    ## Should return list of cross sections of a geometry

    view = GetActiveViewOrCreate('RenderView')
    display = Show(reader, view)
    display.Representation = args['display_representation']
    (xmin,xmax,ymin,ymax,zmin,zmax) = GetActiveSource().GetDataInformation().GetBounds()
    Hide(reader, view)

    slices = []

    for zpos in np.linspace(zmin, zmax, nSlice):
        projection = project(reader, 'Slice', geometry='Plane', origin=[0,0,zpos], normal=[0,0,-1])
        slices.append(projection)

    return slices

def default_origin_normal(reader, origin, normal):
    """
    Input origin and normal parameters for --project should be <float> <string>
    This function takes those values and returns sane vectors to be used in the project() function

    Uses origin as float factor to bounding box limits in the normal direction (unsigned).
    Uses normal direction string to get vectors using direction_handler.

    """
    view = GetActiveViewOrCreate('RenderView')
    display = Show(reader, view)
    (xmin,xmax,ymin,ymax,zmin,zmax) = GetActiveSource().GetDataInformation().GetBounds()
    # print('bbox: ', xmin,xmax,ymin,ymax,zmin,zmax)
    Hide(reader, view)

    new_normal = direction_handler(normal)
    origin_mask = [ xmin + float(origin) * (xmax - xmin), ymin + float(origin) * (ymax - ymin), zmin + float(origin) * (zmax - zmin)]
    new_origin = [abs(x) * y for x,y in zip(new_normal, origin_mask)]

    return new_origin, new_normal

def direction_handler(dir:str):
    """
    Convert from "+x" notation to [1, 0, 0] notation
    """
    dir = dir.strip()
    if len(dir) == 1:
        dir = "+".join(["", dir])

    if dir[1] == "x":
        target = [1, 0, 0]
    elif dir[1] == "y":
        target = [0, 1, 0]
    elif dir[1] == "z":
        target = [0, 0, 1]
    else:
        raise(ValueError)

    if dir[0] == "-":
        target = [-x for x in target]
    elif dir[0] != "+":
        raise(ValueError)

    return target

def view_handler(viewopts:list, zoom:float):
    """
    Set camera view to viewopts ["+x", "-y"]  and zoom
    """
    target = direction_handler(viewopts[0])
    viewup = direction_handler(viewopts[1])

    print("Target:", target)
    print("viewup:", viewup)

    pos = [-x for x in target]

    camera = GetActiveCamera()
    camera.SetFocalPoint(target[0],target[1],target[2])
    camera.SetPosition(pos[0], pos[1], pos[2])
    camera.SetViewUp(viewup[0], viewup[1], viewup[2])
    ResetCamera()

    camera.Dolly(zoom)
    Render()

def main():

    ap = argparse.ArgumentParser()
    ap.add_argument("-cg", "--chromatogram", choices=['Volume', 'Area'], help="Chromatogram: Integrate (0,0,1) surface of given volume or Integrate given area")
    ap.add_argument("-scg", "--shell-chromatograms", type=int, help="Calculate chromatograms in n shell sections of given SURFACE")
    ap.add_argument("--grm-2d", nargs=2, type=int, help="Split into axial and radial sections and integrate scalars for fitting with 2D GRM. args: <ncol> <nrad>")

    ap.add_argument("--integrate", choices=['Volume', 'Area', 'None'], help="Integrate and average the given Volume/Area")


    ap.add_argument("--project", nargs=4, default=['clip', 'Plane', 0.5, "x"], help="Projection. <clip|slice> <Plane|Cylinder..> <origin> <x|y|z>" )
    ap.add_argument("--pipeline", nargs='+', help="Operations to be performed in pipe" )

    ap.add_argument("-st"  , "--shelltype", choices = ['EQUIDISTANT', 'EQUIVOLUME'], default='EQUIDISTANT', help="Shell discretization type")

    ap.add_argument("-sa", "--show-axis", action='store_true', help="Show coordinate axis")
    ap.add_argument("-sb", "--show-scalar-bar", action='store_true', help="Show scalar color bar")
    ap.add_argument("-dr", "--display-representation", default='Surface', choices=['Surface', 'Surface With Edges', 'Points'],  help="Show Surface, Surface With Edges, etc")
    ap.add_argument("-s", "--scalars" , nargs='*' , help="Scalars to consider. (Previously colorvars).")
    ap.add_argument("-z", "--zoom", type=float, default=1, help="Zoom (camera.dolly) value for view")
    ap.add_argument("-v", "--view", nargs=2, default=["+x",  "+y"], help="Set view: target, viewup. Use +x, -z notation.")
    ap.add_argument("-g", "--geometry", nargs=2, type=int, default=[1750, 1300], help="Animation geometry size")
    ap.add_argument("-f", "--filetype", default='pvtu', choices=['xdmf', 'vtu', 'vtk', 'pvtu'], help="filetype: xdmf | vtu | vtk | pvtu")

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

    ## Specifically integrate to get the chromatogram output
    ## with either full volume input -> manually extract surfaces and chromatogram
    ## or surface input, that just needs to be integrated
    integratedData = []
    if args['chromatogram'] == 'Volume':
        GetActiveViewOrCreate('RenderView')

        surfaces = ExtractSurface(Input=reader)
        surfaceNormals = GenerateSurfaceNormals(Input=surfaces)

        threshold = Threshold(Input=surfaceNormals)
        threshold.ThresholdRange = [1.0, 1.0]
        threshold.Scalars = ['POINTS', 'Normals_Z']

        integratedData = integrate(threshold, scalars, normalize='Area', timeArray=timeArray)

        for scalar in scalars:
            csvWriter("".join([scalar, '.csv']), timeArray, np.array(integratedData).T[list(scalars).index(scalar)])

    elif args['chromatogram'] == 'Area':
        integratedData = integrate(reader, scalars, normalize='Area', timeArray=reader.TimestepValues)

        for scalar in reader.PointArrayStatus:
            csvWriter("".join([scalar, '.csv']), reader.TimestepValues, np.array(integratedData).T[list(reader.PointArrayStatus).index(scalar)])

    elif args['shell_chromatograms']:
        print("Running shell_chromatograms on provided SURFACE!")
        nRegions = int(args['shell_chromatograms'])

        view = GetActiveViewOrCreate('RenderView')
        ## Calc bounding box. Requires show
        display = Show(reader, view)
        (xmin,xmax,ymin,ymax,zmin,zmax) = GetActiveSource().GetDataInformation().GetBounds()
        Hide(reader, view)

        nShells = nRegions + 1 #Including r = 0
        rShells = []

        R = (xmax - xmin + ymax - ymin)/4
        print("R:", R)

        # shellType = 'EQUIDISTANT'
        # shellType = 'EQUIVOLUME'
        if shellType == 'EQUIVOLUME':
            for n in range(nShells):
                rShells.append(R * sqrt(n/nRegions))
        elif shellType == 'EQUIDISTANT':
            for n in range(nShells):
                rShells.append(R * (n/nRegions))

        print("rShells:", rShells)

        radAvg = []
        integrated_over_time = [ [] for region in range(nRegions) ]
        for timestep in range(nts):

            timeKeeper.Time = timestep
            reader.UpdatePipeline(reader.TimestepValues[timestep])

            print("its:", timestep)

            for radIn, radOut in zip(rShells[:-1], rShells[1:]):

                index = rShells.index(radIn)

                radAvg.append( (radIn + radOut) / 2 )

                shell_area = np.pi * (radOut**2 - radIn**2)

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

                ## Since we handle time outside the integrate function, the
                ## only entry in the outer list is the list of integrated scalars
                ## at the given time
                integratedData = integrate(clipInner, scalars, normalize='Area')
                integrated_over_time[index].extend(integratedData)

                Delete(clipInner)
                Delete(clipOuter)

        for region in range(nRegions):
            for scalar in scalars:
                csvWriter("shell_{i}_{s}.cg".format(i=region, s=scalar), timeArray, np.array(integrated_over_time[region]).T[list(scalars).index(scalar)])

    elif args['grm_2d']:
        GRM2D(reader, args)
    elif args['integrate']:
        integrated_over_time = integrate(reader, args['scalars'], normalize=args['integrate'], timeArray=timeArray)
        print(integrated_over_time)
        for scalar in scalars:
            csvWriter("{s}.integrated.{n}.csv".format(s=scalar, n=args['integrate']), timeArray, np.array(integrated_over_time).T[list(scalars).index(scalar)])


    ## NOTE: Pipeline operations below [EXPERIMENTAL]
    ## The idea is to provide the sequence of operations on the commandline
    ##  and execute it here
    supported_operations = {
        'project': project,
        'screenshot': screenshot
    }

    object = reader
    args['pipeline'] = args['pipeline'] or []

    for operation in args['pipeline']:
        object = supported_operations[operation](object, args)

if __name__ == "__main__":
    main()
