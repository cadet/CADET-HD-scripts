#!/bin/env -S pvpython --force-offscreen-rendering

""" The successor script to vis.py

    - More modular
    - Built for generic XNS

"""

from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()
from vtkmodules.numpy_interface import dataset_adapter as dsa
import vtk.util.numpy_support as ns #type:ignore
from math import sqrt
import numpy as np

import argparse
import os
import csv

def csvWriter(filename, x, y):
    with open(filename, 'w') as f:
        writer = csv.writer(f)
        writer.writerows(zip(x, y))

def setCameraOrientation(zoom):
    camera = GetActiveCamera()
    camera.SetFocalPoint(-1,0,0)
    camera.SetPosition(1,0,0)
    camera.SetViewUp(0,-1,0)
    ResetCamera()

    camera.Dolly(zoom)
    Render()

def integrate(object, vars, normalize=None, timeArray=[]):
    ## normalize= "Volume" or "Area" or None
    choices = ['Volume', 'Area']

    timeKeeper = GetTimeKeeper()
    nts = len(timeArray) or 1

    # GetActiveViewOrCreate('RenderView')


    integrated_over_time = []
    for timestep in range(nts):
        try:
            timeKeeper.Time = timestep
            object.UpdatePipeline(timeArray[timestep])
        except IndexError:
            pass

        integrated = IntegrateVariables(Input=object)
        intdata = servermanager.Fetch(integrated)
        intdata = dsa.WrapDataObject(intdata)

        volume=1
        if normalize in choices:
            if normalize in intdata.CellData.keys():
                volume = intdata.CellData[normalize][0]
            else:
                print("".join(["Cannot normalize by ", normalize, ". No such CellData!"]))

        print("{key}: {value}".format(key=normalize, value=volume))

        integrated_scalars = []
        for var in vars:
            value = intdata.PointData[var]
            value = ns.vtk_to_numpy(value)
            integrated_scalars.append(value[0]/volume)  ## Average of c, instead of integ(c.dV)

    # return integrated_scalars

        integrated_over_time.append(integrated_scalars)
    return integrated_over_time

def project(inputView, projectionType, geometry='Plane', origin=None, normal=[1,0,0]):

    projectionView = GetActiveViewOrCreate('RenderView')
    display = Show(inputView, projectionView)

    (xmin,xmax,ymin,ymax,zmin,zmax) = GetActiveSource().GetDataInformation().GetBounds()
    Hide(inputView, projectionView)

    center = [ (xmax+xmin)/2, (ymax+ymin)/2, (zmax+zmin)/2,]

    if not origin:
        origin=center

    projection = None
    if projectionType == 'Clip':
        projection = Clip(Input=inputView)
        projection.ClipType = geometry
        projection.HyperTreeGridSlicer = geometry
        projection.ClipType.Origin = origin
        projection.ClipType.Normal = normal
        Hide3DWidgets(proxy=projection.ClipType)
    elif projectionType == 'Slice':
        projection = Slice(Input=inputView)
        projection.SliceType = geometry
        projection.HyperTreeGridSlicer = geometry
        projection.SliceType.Origin = origin
        projection.SliceType.Normal = normal
        Hide3DWidgets(proxy=projection.SliceType)
    else:
        projection = inputView

    projection.UpdatePipeline()

    return projection

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

# def get_radial_shells(reader):

def main():

    ap = argparse.ArgumentParser()
    ap.add_argument("-cg", "--chromatogram", help="Chromatogram: Integrate (0,0,1) surface of given volume, Integrate Area")
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
        sys.exit(-1)

    timeKeeper = GetTimeKeeper()
    timeArray = reader.TimestepValues

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

        integratedData = integrate(threshold, reader.PointArrayStatus , normalize='Area', timeArray=timeArray)

        for scalar in reader.PointArrayStatus:
            csvWriter("".join([scalar, '.csv']), reader.TimestepValues, np.array(integratedData).T[list(reader.PointArrayStatus).index(scalar)])

    elif args['chromatogram'] == 'Area':
        integratedData = integrate(reader, reader.PointArrayStatus, normalize='Area', timeArray=reader.TimestepValues)
        print(integratedData)

        for scalar in reader.PointArrayStatus:
            csvWriter("".join([scalar, '.csv']), reader.TimestepValues, np.array(integratedData).T[list(reader.PointArrayStatus).index(scalar)])

if __name__ == "__main__":
    main()
