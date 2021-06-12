#!/bin/env -S pvpython --force-offscreen-rendering

""" The successor script to vis.py

    - More modular
    - Built for generic XNS

"""

def setCameraOrientation(zoom):
    camera = GetActiveCamera()
    camera.SetFocalPoint(-1,0,0)
    camera.SetPosition(1,0,0)
    camera.SetViewUp(0,-1,0)
    ResetCamera()

    camera.Dolly(zoom)
    Render()

def integrate(object, nvars, normalize=None):
    ## normalize= "Volume" or "Area" or None

    integrated = IntegrateVariables(Input=object)
    intdata = servermanager.Fetch(integrated)
    intdata = dsa.WrapDataObject(intdata)

    volume=1
    if normalize:
        volume = intdata.CellData[normalize][0]

    values = []
    for i in range(nvars):
        value = intdata.PointData[i]
        value = ns.vtk_to_numpy(value)
        values.append(value[0]/volume)  ## Average of c, instead of integ(c.dV)

    return values

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

