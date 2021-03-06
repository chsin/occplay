import sys
from OCC.Display.SimpleGui import *
from OCC.gp import *
from OCC.BRepPrimAPI import *
from OCC.BRepAlgoAPI import *
from OCC.BRepFilletAPI import *
from OCC.TopExp import *
from OCC.TopoDS import *
from OCC.TopAbs import *
from OCC import StlAPI
from OCC import STEPControl
from OCC.Utils import Topology
from OCC.AIS import * 
import numpy as np
from scipy import linalg

reScale = 3.93701  
cubeSide = 300.0 * 1.5 / 2.0 / reScale
sphereRadius = 100.0 * 1.5 / 2.0 / reScale
tanAngle = 1.5
apexShift = sphereRadius * (np.sqrt(1.0 - (4.0/13.0)**2) / tanAngle - 4.0/13.0)
cylinderHeight = 12.8 * 2 / reScale
cylinderDiameter = 25.6 / reScale


def toSTL(shape, filename, ascii=False, deflection=0.001):
    stl_writer = StlAPI.StlAPI_Writer()
    stl_writer.SetASCIIMode(ascii)
    stl_writer.SetDeflection(deflection)
    stl_writer.Write(shape, filename)


def toSTEP(shape, filename, verbose=False, tolerance=0.001):
    '''
    Writes STL output of the solid

    :param filename: Path of the file to write JSON to
    :type filename: str
    :param verbose: Choose if you want to see the STEP stats
    :type verbose: bool
    :param tolerance: Provides control over quality of exported STEP.
    :type tolernace: float
    '''
    class Verboseness():
      def __enter__(self):
        if not verbose:
          import sys, os
          sys.stdout.flush()
          self.newstdout = os.dup(1)
          self.devnull = os.open('/dev/null', os.O_WRONLY)
          os.dup2(self.devnull, 1)
          os.close(self.devnull)
          sys.stdout = os.fdopen(self.newstdout, 'w')
      def __exit__(self, type, value, traceback):
        if not verbose:
          os.dup2(self.newstdout, 1)

    with Verboseness():
      stepWriter = STEPControl.STEPControl_Writer()
      stepWriter.SetTolerance(tolerance)
      if shape:
        status = stepWriter.Transfer(shape, 
                                     STEPControl.STEPControl_AsIs)
        if status:
          stepWriter.Write(filename)


def make_point(coordinates):
    """ Creates a point given its coordinates as numpy.array. """
    # Components in floats
    x = float(coordinates[0])
    y = float(coordinates[1])
    z = float(coordinates[2])
    return gp_Pnt(x, y, z)

def make_unit_vector(direction):
    """  Creates a unit vector given its direction as numpy.array. """
    # Normalize the direction
    direction /= linalg.norm(direction)
    # Components in floats
    x = float(direction[0])
    y = float(direction[1])
    z = float(direction[2])
    return gp_Dir(x, y, z)


def make_axis(apex, axis):
    """  Creates an axis --- a tethered unit vector given its start point
    and direction as numpy.arrays. """
    origin = make_point(apex)
    direction = make_unit_vector(axis)
    return gp_Ax2(origin, direction)


def make_sphere(center, radius):
    '''
    Creates a sphere given center and radius.
    @param center: center of a sphere as a numpy array
    @type center: numpy.ndarray(3,)
    @param radius: radius of the sphere
    @type radius: float
    '''
    point = make_point(center)
    return BRepPrimAPI_MakeSphere(point, radius)


def make_cylinder(apex, axis, height, radius):
    """
    Creates a cylinder.
    @param apex: start point of cylinder
    @type apex: numpy.ndarray(3,)
    @param axis: direction vector of the cylinder's main axis
    @type directionvector: numpy.ndarray(3,)
    @param height: height of cylinder
    @type height: float
    @param radius: radius of cylinder
    @type radius: float
    """
    cylinderAxis = make_axis(apex, axis)
    cylinder = BRepPrimAPI_MakeCylinder(cylinderAxis, radius, height)
    return cylinder


def make_cone(apex, axis, height, radius1, radius2):
    """
    Creates a cone or a frustum.
    param apex: start point of cone
    @type apex: numpy.ndarray(3,)
    @param axis: direction vector of the cone's main axis
    @type directionvector: numpy.ndarray(3,)
    @param height: cone height
    @type height: float
    @param radius1: radius at the cone bottom
    @type radius1: float
    @param radius1: radius at the cone top
    @type radius1: float
    """
    coneAxis = make_axis(apex, axis)
    cone = BRepPrimAPI_MakeCone(coneAxis, radius1, radius2, height)
    return cone


def make_cube( side ):
    '''
    Creates a cube given its side.
    @param side: side of the cube
    @type side: float
    '''
    a = float(side / 2.0);
    ones = np.ones(3, dtype=float)
    return BRepPrimAPI_MakeBox(make_point(-a * ones), make_point( a * ones))


def layer_shape(side, v, sign = 1):
    a = float(side / 2.0) * sign
    array1 = -a * np.ones(3, dtype=float)
    point1 = make_point(-a * np.ones(3, dtype=float))
    array2 =  a *(np.ones(3, dtype=float) - 4./3 * v)
    point2 = make_point(a *(np.ones(3, dtype=float) - 4./3 * v))
    return BRepPrimAPI_MakeBox(point1, point2).Shape()


def sphere_shape():
    origin = np.zeros(3, dtype=float)
    sphere = make_sphere( origin, sphereRadius )
    return sphere.Shape()


def cylinder_shape(axis):
    apex = (sphereRadius - cylinderHeight) * axis / linalg.norm(axis)
    radius = cylinderDiameter / 2.0
    return make_cylinder(apex, axis, cylinderHeight * 2.0, radius).Shape()


def cone_shape(axis):
    # create cone with this axis
    apex = apexShift * axis
    height = cubeSide
    radius = tanAngle * height
    return make_cone(apex, -axis, height, 0.0, radius).Shape()


def cube_shape():
    return make_cube(cubeSide).Shape()


def fillet_all(shape):
    fillet = BRepFilletAPI_MakeFillet(shape)
    for edge in Topology.Topo(shape).edges():
        fillet.Add(.5, edge)
    return fillet.Shape()


def draw_sphere(event = None):
    display.DisplayColoredShape( sphere_shape() , 'RED' )


def draw_cube(event=None):
    display.DisplayColoredShape( cube_shape() , 'WHITE' )


def cut_bottom_cone():
    coneShape = cone_shape(np.array((-1, 0, 0), dtype=np.float))
    cone_minus_sphere = BRepAlgoAPI_Cut(coneShape, sphere_shape()).Shape()
    return BRepAlgoAPI_Common(cone_minus_sphere, cube_shape()).Shape()


def cut_top_cone():
    coneShape = cone_shape(np.array((1, 0, 0), dtype=np.float))
    cone_minus_sphere = BRepAlgoAPI_Cut(coneShape, sphere_shape()).Shape()
    return BRepAlgoAPI_Common(cone_minus_sphere, cube_shape()).Shape()


def cone_layer_shape(direction, sign = 1.0):
    coneShape = cone_shape(direction * sign)
    layer = layer_shape(cubeSide, direction, sign)
    step1 =  BRepAlgoAPI_Common(coneShape, cube_shape()).Shape()
    step2 =  BRepAlgoAPI_Fuse(step1, layer).Shape()
    return BRepAlgoAPI_Cut(step2, sphere_shape()).Shape()


def x_cone_layer(sign = 1.0):
    return cone_layer_shape(np.array((1, 0, 0), dtype=np.float), sign)


def y_cone_layer(sign = 1.0):
    return cone_layer_shape( np.array((0, 1, 0), dtype=np.float), sign)


def z_cone_layer(sign = 1.0):
    return cone_layer_shape(np.array((0, 0, 1), dtype=np.float), sign)


def face_piece():
    step1 = BRepAlgoAPI_Cut(y_cone_layer(-1.0), z_cone_layer(-1.0)).Shape()
    step2 = BRepAlgoAPI_Cut(step1, z_cone_layer(1.0)).Shape()
    step3 = BRepAlgoAPI_Cut(step2, x_cone_layer(1.0)).Shape()
    step4 = BRepAlgoAPI_Cut(step3, x_cone_layer(-1.0)).Shape()
    v = np.array((0, 1, 0), dtype=np.float)
    return BRepAlgoAPI_Cut(step4, cylinder_shape(v)).Shape()


def edge_piece():
    step1 = BRepAlgoAPI_Common(y_cone_layer(-1.0), z_cone_layer(-1.0)).Shape()
    step2 = BRepAlgoAPI_Cut(step1, x_cone_layer(1.0)).Shape()
    step3 = BRepAlgoAPI_Cut(step2, x_cone_layer(-1.0)).Shape()
    v = np.array((0, 1, 1), dtype=np.float)
    return BRepAlgoAPI_Cut(step3, cylinder_shape(v)).Shape()


def corner_piece():
    step1 = BRepAlgoAPI_Common(y_cone_layer(-1.0), z_cone_layer(-1.0)).Shape()
    step2 = BRepAlgoAPI_Common(step1, x_cone_layer(1.0)).Shape()
    v = np.array((-1, 1, 1), dtype=np.float)
    return BRepAlgoAPI_Cut(step2, cylinder_shape(v)).Shape()

def holder():
    #variables used
    stock2D = 300.0 * 1.0 / reScale #3"
    cubeSidePlus = cubeSide*1.004  #2.259"
    stockHalved = float(stock2D / 2.0) # 1.5"
    csp3rd = float(cubeSidePlus / 3.0) # 0.753"

    #make stock: 3x3x.753
    stockP1 = make_point(np.array((-stockHalved, -stockHalved, 0), dtype=np.float))
    stockP2 = make_point(np.array((stockHalved, stockHalved, csp3rd), dtype=np.float))
    stock_box = BRepPrimAPI_MakeBox(stockP1, stockP2).Shape()

    #cut where rubiks layer will sit:
    cb_XY = cubeSidePlus / 2.0
    cbZ_low = cubeSide / 2.0 / 3.0
    cbZ_hi = cubeSide / 2.0
    cbP1 = make_point(np.array((-cb_XY, -cb_XY, cbZ_low), dtype=np.float))
    cbP2 = make_point(np.array((cb_XY, cb_XY, cbZ_hi), dtype=np.float))
    rubiks_layer_box = BRepPrimAPI_MakeBox(cbP1, cbP2).Shape()

    #corners
    dir_pnt = np.array((0, 0, 1), dtype=np.float)
    r = 100.0 * .25 / 2.0 / reScale #.125
    h = cbZ_hi - cbZ_low
    xy_list = [(cb_XY,cb_XY), (cb_XY, -cb_XY), (-cb_XY,cb_XY), (-cb_XY,-cb_XY)]
    for xy in xy_list:
        cpt = np.array((xy[0], xy[1], cbZ_low), dtype=np.float)
        cc = make_cylinder(cpt, dir_pnt, h, r).Shape()
        rubiks_layer_box = BRepAlgoAPI_Fuse(rubiks_layer_box, cc).Shape()

    # remove insides
    cutRLayer = BRepAlgoAPI_Cut(stock_box, rubiks_layer_box).Shape()

    # cut cylinder for cubes' conical protrusion
    sphereRadius = 100.0 * 1.5 / 2.0 / reScale #d=1.5,r=.75
    center_pnt = np.array((0, 0, sphereRadius*.25), dtype=np.float)
    cylinder = make_cylinder(center_pnt, dir_pnt, csp3rd, sphereRadius*1.3).Shape()
    holder = BRepAlgoAPI_Cut(cutRLayer, cylinder).Shape()

    # add sphere shape
    bottomP1 = make_point(np.array((-stockHalved, -stockHalved, -stockHalved), dtype=np.float))
    bottomP2 = make_point(np.array((stockHalved, stockHalved, 0), dtype=np.float))
    box_bottom = BRepPrimAPI_MakeBox(bottomP1, bottomP2).Shape()
    cutted_sphere = BRepAlgoAPI_Cut(sphere_shape(), box_bottom).Shape()
    addSphere = BRepAlgoAPI_Fuse(cutted_sphere, holder).Shape()

    # displays shape
    ais_shape = AIS_Shape(addSphere).GetHandle()
    ais_context = display.GetContext().GetObject()
    ais_context.SetTransparency(ais_shape, 0.5, True)
    ais_context.Display(ais_shape)

#### Draw functions
def draw_cut_top_cone(event=None):
    display.DisplayColoredShape(cut_top_cone(), 'BLACK')


def draw_cut_bottom_cone(event=None):
    display.DisplayColoredShape(cut_bottom_cone(), 'BLACK')


def draw_cylinder(event=None):
    v = np.array((-1, 1, 1), dtype=np.float)
    display.DisplayColoredShape(cylinder_shape(v), 'WHITE')

def draw_white_layer(event=None):
    display.DisplayColoredShape(x_cone_layer(1.0), 'BLUE')

def draw_yellow_layer(event=None):
    display.DisplayColoredShape(x_cone_layer(-1.0), 'BLUE')


def draw_blue_layer(event=None):
    display.DisplayColoredShape( y_cone_layer(1.0), 'BLUE')


def draw_green_layer(event=None):
    display.DisplayColoredShape(y_cone_layer(-1.0), 'BLUE')


def draw_red_layer(event=None):
    display.DisplayColoredShape(z_cone_layer(1.0), 'BLUE')


def draw_orange_layer(event=None):
    display.DisplayColoredShape(z_cone_layer(-1.0), 'BLUE')


def draw_face_piece(event=None):
    display.DisplayColoredShape(fillet_all(face_piece()), 'BLACK')


def draw_edge_piece(event=None):
    display.DisplayColoredShape(fillet_all(edge_piece()), 'GREEN')

def draw_corner_piece(event=None):
    print corner_piece()
    display.DisplayColoredShape(fillet_all(corner_piece()), 'CYAN')

def erase_all(event=None):
    display.EraseAll()


if __name__ == '__main__':
#    toSTEP(corner_piece(), 'corner.step')
#    toSTEP(edge_piece(), 'edge.step')
#    toSTEP(face_piece(), 'face.step')

#    toSTL(corner_piece(), 'corner.stl')
#    toSTL(edge_piece(), 'edge.stl')
#    toSTL(face_piece(), 'face.stl')

    display, start_display, add_menu, add_function_to_menu = \
        init_display()
    add_menu('Draw')
    add_function_to_menu('Draw', draw_sphere)
    add_function_to_menu('Draw', draw_cylinder)
    add_function_to_menu('Draw', draw_cube)
    add_function_to_menu('Draw', draw_cut_bottom_cone)
    add_function_to_menu('Draw', draw_cut_top_cone)
    add_menu('Layers')
    add_function_to_menu('Layers', draw_white_layer)
    add_function_to_menu('Layers', draw_yellow_layer)
    add_function_to_menu('Layers', draw_green_layer)
    add_function_to_menu('Layers', draw_red_layer)
    add_function_to_menu('Layers', draw_blue_layer)
    add_function_to_menu('Layers', draw_orange_layer)
    add_menu('Cubes')
    add_function_to_menu('Cubes', draw_face_piece)
    add_function_to_menu('Cubes', draw_edge_piece)
    add_function_to_menu('Cubes', draw_corner_piece)
    add_menu('Holder')
    add_function_to_menu('Holder', holder)
    add_menu('Erase')
    add_function_to_menu('Erase', erase_all)
    start_display()
