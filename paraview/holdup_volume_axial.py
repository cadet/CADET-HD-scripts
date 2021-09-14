#!/usr/bin/env pvpython

"""
Script to calculate holdup volumes at different axial lengths of a given column.

Assumes input to be concentrations in only the interstitial region.
"""

# import numpy as np
import math
import argparse
import os
from paraview.simple import *
from vtkmodules.numpy_interface import dataset_adapter as dsa
import vtk.util.numpy_support as ns #type:ignore

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

def ana_holdup_vol(v_i, v_b):
    """
    Calculate analytical holdup volume
    """
    eps  = 0.75
    qmax = 4.88
    ka   = 1.144
    kd   = 2.0e-3
    cin  = 7.14e-3
    qinf = qmax * ka / (ka*cin + kd)
    holdup_real = (v_i)  + (eps * v_b)  +  ( (1-eps) * v_b * qinf )
    return holdup_real

def trapz(y, *args, **kwargs):
    """
    trapezoidal integration of a curve. Originally created to avoid using numpy with postchrom.py on JURECA.
    """
    x = kwargs.get('x', range(len(y)))
    sum = 0.0
    for i in range(len(x) - 1):
     sum = sum + 0.5 * (y[i] + y[i+1]) * (x[i+1] - x[i])

    return sum

def num_holdup_vol(t, c, R, u, cin):
    """
    calculate numerical holdup volume from a chromatogram (t, c)
    """
    cn = [ (1 - elem / cin) for elem in c]
    holdup_num = trapz(cn, x=t) * math.pi * R**2 * u
    return holdup_num

def project(object, args):

    projectionType = args['project'][0]
    geometry = args['project'][1]
    origin = args['project'][2]
    normal = args['project'][3]

    origin, normal = default_origin_normal(object, origin, normal)

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

        if not intdata:
            raise(AssertionError)

        volume=1
        if normalize in choices:
            if normalize in intdata.CellData.keys():
                volume = intdata.CellData[normalize][0]
            else:
                print("".join(["Cannot normalize by ", normalize, ". No such CellData!"]))

        # print("{key}: {value}".format(key=normalize, value=volume))

        integrated_scalars = []
        for var in vars:
            value = intdata.PointData[var]
            value = ns.vtk_to_numpy(value)
            integrated_scalars.append(value[0]/volume)  ## Average of c, instead of integ(c.dV)

        integrated_over_time.append(integrated_scalars)

        Delete(integrated)

    return integrated_over_time

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project", nargs=4, default=['clip', 'Plane', 0.5, "x"], help="Projection. <clip|slice> <Plane|Cylinder..> <origin> <x|y|z>" )

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

    #######

    view = GetActiveViewOrCreate('RenderView')
    ## Calc bounding box. Requires show
    display = Show(reader, view)
    (xmin,xmax,ymin,ymax,zmin,zmax) = GetActiveSource().GetDataInformation().GetBounds()
    Hide(reader, view)

    R_cyl = ((xmax - xmin) + (ymax - ymin))/ 4
    length_full = zmax - zmin

    cut_fraction = 0.1

    clipped = project(reader, { 'project': ['clip', 'Plane', cut_fraction, '+z'] })
    sliced  = project(reader, { 'project': ['slice', 'Plane', cut_fraction, '+z'] })

    integrated = IntegrateVariables(Input=clipped)
    intdata = servermanager.Fetch(integrated)
    intdata = dsa.WrapDataObject(intdata)
    volume_int = intdata.CellData['Volume'][0]

    volume_cyl = math.pi * R_cyl**2 * length_full * cut_fraction

    volume_beads = volume_cyl - volume_int

    HV_ana = ana_holdup_vol(volume_int, volume_beads)

    integratedData = integrate(sliced, scalars, normalize='Area', timeArray=timeArray)
    integratedData = [ y for x in integratedData for y in x ]

    print(integratedData)
    print(timeArray)

    HV_num = num_holdup_vol(timeArray, integratedData, R_cyl, 2.09e-4, 7.14e-3)

    print("HV_ana =", HV_ana)
    print("HV_num =", HV_num)

if __name__ == "__main__":
    main()
