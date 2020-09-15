#!/bin/env pvpython

## TODO: snapshot time range
## TODO: snapshot time stride
## DONE: distribute to parallel

import argparse
from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()
from vtkmodules.numpy_interface import dataset_adapter as dsa
import vtk.util.numpy_support as ns
import numpy as np
from matplotlib import pyplot as plt

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
    colorVar = kwargs.get('colorVar', reader.PointArrayStatus[0])

    animationScene1 = GetAnimationScene()
    timeKeeper1 = GetTimeKeeper()
    animationScene1.UpdateAnimationUsingDataTimeSteps()

    ## Use last timestep as reference for creating color map
    animationScene1.AnimationTime = reader.TimestepValues[-1]
    timeKeeper1.Time = reader.TimestepValues[-1]

    renderView1 = GetActiveViewOrCreate('RenderView')
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


    projectionDisplay = Show(projection, renderView1)
    projectionDisplay.Representation = 'Surface'

    ColorBy(projectionDisplay, ('POINTS', colorVar))

    wLUT = GetColorTransferFunction(colorVar)
    wPWF = GetOpacityTransferFunction(colorVar)

    projectionDisplay.UpdatePipeline()
    projectionDisplay.RescaleTransferFunctionToDataRange(True,False)
    projectionDisplay.SetScalarBarVisibility(renderView1, False)

    renderView1.Update()
    renderView1.ResetCamera()

    renderView1.CameraPosition = [0.0005945160428284565, 1.959300980464672e-13, -1.3552527156068805e-20]
    renderView1.CameraFocalPoint = [-2.549999918319142e-05, 1.959300980464672e-13, -1.3552527156068805e-20]
    renderView1.CameraViewUp = [0.0, 0.0, 1.0]
    renderView1.CameraParallelScale = 0.00016047195994169908
    renderView1.ResetCamera()

    # SaveScreenshot('FLOW_'+colorVar+'.png', renderView1, ImageResolution=[1750, 1300], TransparentBackground=1)

    SaveAnimation(colorVar + '.png', renderView1, ImageResolution=[1750, 1300], TransparentBackground=1)


def main():

    ap = argparse.ArgumentParser()
    ap.add_argument("-c", "--colorVar", required=False,
            help="color map variable")

    ap.add_argument("-p", "--projectionType", required=False, default='clip',
            help="projection type: clip | slice")
    ap.add_argument("-s", "--snapshot", required=False, action='store_true',
            help="run snapshotter")
    ap.add_argument("-d", "--distribute", required=False,
            help="Apply d3 filter and save data")
    ap.add_argument("-m", "--mass-defect", required=False, action='store_true',
            help="Integrate and find mass_defect curve")

    ap.add_argument("-r", "--reader", required=True,
            help="reader: xdmf | vtu | vtk | pvtu")
    ap.add_argument("-w", "--writer", required=False,
            help="reader: xdmf | vtu | vtk | pvtu")

    ap.add_argument("FILES", nargs='*',
            help="files..")
    args = vars(ap.parse_args())

    for key in args:
        print(key + ': ', args[key])

    if not args['FILES']:
        print("No Input Files Given!")
        sys.exit(0)

    reader=None
    if args['reader'] == 'xdmf':
        reader = XDMFReader(FileNames=args['FILES'])
    elif args['reader'] == 'vtu':
        reader = XMLUnstructuredGridReader(FileName=args['FILES'])
    elif args['reader'] == 'pvtu':
        reader = XMLPartitionedUnstructuredGridReader(FileName=args['FILES'])
    elif args['reader'] == 'vtk':
        reader = LegacyVTKReader(FileName=args['FILES'])
    else:
        print("Reader Unspecified!")
        sys.exit(-1)

    if args['distribute']:
        d3 = D3(Input=reader)
        SaveData(args['distribute'] + '.pvtu', proxy=d3,
                Writealltimestepsasfileseries=1)
        sys.exit(0)

    if args['snapshot']:
        snapshot(reader,
                projectionType = args['projectionType'],
                colorVar = args['colorVar'])

    if args['mass_defect']:
        mass_defect(reader);

    if args['writer']:
        print("ENTERED WRITER")
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
