import sys
from OCC.Display.SimpleGui import *
from OCC.gp import *
from OCC.BRepPrimAPI import *
from OCC.BRepAlgoAPI import *
from OCC.StlAPI import *
from OCC.STEPControl import *
import numpy as np
from scipy import linalg

cubeSide = 300.0
sphereDiameter = 200.0

def toSTL(shape, filename, ascii=False, deflection=0.001):
    stl_writer = StlAPI_Writer()
    stl_writer.SetASCIIMode(ascii)
    stl_writer.SetDeflection(deflection)
    stl_writer.Write(shape, filename)

def toSTEP(shape, filename, verbose=False, tolerance=0.001):
        stepWriter = STEPControl_Writer()
        stepWriter.SetTolerance(tolerance)
        if shape:
            status = stepWriter.Transfer(shape, STEPControl_AsIs )
        if status:
            stepWriter.Write(filename)
        if verbose:
            stepWriter.PrintStatsTransfer()

def make_point(coordinates):
    '''
    Creates a point given its coordinates as numpy.array.
    @param coordinates: coordinates as a numpy array
    @type coordinates: numpy.ndarray(3,)
    '''
    # Components in floats
    x = float(coordinates[0])
    y = float(coordinates[1])
    z = float(coordinates[2])
    return gp_Pnt(x, y, z)

def make_unit_vector(direction):
    '''
    Creates a unit vector given its direction as numpy.array.
    @param direction: direction as a numpy array
    @type direction: numpy.ndarray(3,)
    '''
    # Normalize the direction
    direction /= linalg.norm(direction)
    # Components in floats
    x = float(direction[0])
    y = float(direction[1])
    z = float(direction[2])
    return gp_Dir(x, y, z)

def make_axis(apex, axis):
    '''
    Creates an axis --- a tethered unit vector given its start point
    and direction as numpy.arrays.
    @param apex: start poing coordinates as a numpy array
    @type apex: numpy.ndarray(3,)
    @param axis: direction as a numpy array
    @type axis: numpy.ndarray(3,)
    '''
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

def make_cone(apex, axis, height, radius1, radius2):
    """
    Creates a cone or a frustum.
    @param vector: vector at the beginning
    @type vector: numpy.ndarray(3,)
    @param directionvector: direction vector of the cone maina axis
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
    sphere = make_sphere( origin, sphereDiameter / 2.0 )
    return sphere.Shape()
   
def cone_shape(v):
    # create cone with axis in v directon
    apex = np.array(20 * v, dtype=float)
    axis = np.array(-100 * v, dtype=float)
    height = 300.0
    radius1 = 0.0
    radius2 = 600.0
    return make_cone(apex, axis, height, radius1, radius2).Shape()

def cube_shape():
    return make_cube(cubeSide).Shape()

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
    return BRepAlgoAPI_Cut(step3, x_cone_layer(-1.0)).Shape()

def edge_piece():
    step1 = BRepAlgoAPI_Common(y_cone_layer(-1.0), z_cone_layer(-1.0)).Shape()
    step2 = BRepAlgoAPI_Cut(step1, x_cone_layer(1.0)).Shape()
    return BRepAlgoAPI_Cut(step2, x_cone_layer(-1.0)).Shape()

def corner_piece():
    step = BRepAlgoAPI_Common(y_cone_layer(-1.0), z_cone_layer(-1.0)).Shape()
    return BRepAlgoAPI_Common(step, x_cone_layer(1.0)).Shape()

def draw_cut_top_cone(event=None):
    display.DisplayColoredShape(cut_top_cone(), 'BLACK')

def draw_cut_bottom_cone(event=None):
    display.DisplayColoredShape(cut_bottom_cone(), 'BLACK')

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
    display.DisplayColoredShape(face_piece(), 'CYAN')

def draw_edge_piece(event=None):
    display.DisplayColoredShape(edge_piece(), 'GREEN')

def draw_corner_piece(event=None):
    display.DisplayColoredShape(corner_piece(), 'BLACK')

def erase_all(event=None):
    display.EraseAll()

if __name__ == '__main__':
    toSTEP(corner_piece(), 'corner.step')
    toSTEP(edge_piece(), 'edge.step')
    toSTEP(face_piece(), 'face.step')

    toSTL(corner_piece(), 'corner.stl')
    toSTL(edge_piece(), 'edge.stl')
    toSTL(face_piece(), 'face.stl')

    display, start_display, add_menu, add_function_to_menu = \
        init_display()
    add_menu('Draw')
    add_function_to_menu('Draw', draw_sphere)
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
    add_menu('Erase')
    add_function_to_menu('Erase', erase_all)
    start_display()
