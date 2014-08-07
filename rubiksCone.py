import sys
from OCC.Display.SimpleGui import *
from OCC.gp import *
from OCC.BRepPrimAPI import *
from OCC.BRepAlgoAPI import *
from OCC.StlAPI import *
import numpy as np
from scipy import linalg

def toSTL(thingus, filename, ascii=False, deflection=0.001):
    stl_writer = StlAPI_Writer()
    stl_writer.SetASCIIMode(ascii)
    stl_writer.SetDeflection(deflection)
    stl_writer.Write(thingus, filename)

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

def sphere_from_center_and_radius( center,
                                   radius ):
    '''
    Creates a sphere given center and radius.
    @param center: center of a sphere as a numpy array
    @type center: numpy.ndarray(3,)
    @param radius: radius of the sphere
    @type radius: float
    '''
    point = make_point(center)
    return BRepPrimAPI_MakeSphere( point, radius )

def cone_from_point_height_directionvector_and_two_radii( apex,
                                                          axis,
                                                          height,
                                                          radius1,
                                                          radius2 ):
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
    cone = BRepPrimAPI_MakeCone( coneAxis,
                                 radius1,
                                 radius2,
                                 height )
    return cone

def centered_cube_from_side( side ):
    '''
    Creates a cube given its side.
    @param side: side of the cube
    @type side: float
    '''
    a = float(side / 2.0);
    point1 = make_point(-a * np.ones(3, dtype=float))
    point2 = make_point( a * np.ones(3, dtype=float))
    return BRepPrimAPI_MakeBox(point1, point2)

def top_layer( side ):
    a = float(side / 2.0)
    array1 = -a * np.ones(3, dtype=float)
    point1 = make_point(array1)
    array2 = np.array((-a/3.0 , a, a), dtype=float)
    point2 = make_point(array2) 
    return BRepPrimAPI_MakeBox(point1, point2)

def bottom_layer( side ):
    a = float(side / 2.0)
    array1 = a * np.ones(3, dtype=float)
    point1 = make_point(array1)
    array2 = np.array((a/3.0 , -a, -a), dtype=float)
    point2 = make_point(array2) 
    return BRepPrimAPI_MakeBox(point1, point2)

def left_layer( side ):
    a = float(side / 2.0)
    array1 = -a * np.ones(3, dtype=float)
    point1 = make_point(array1)
    array2 = np.array((a, a, -a / 3.0), dtype=float)
    point2 = make_point(array2) 
    return BRepPrimAPI_MakeBox(point1, point2)

def right_layer( side ):
    a = float(side / 2.0)
    array1 = a * np.ones(3, dtype=float)
    point1 = make_point(array1)
    array2 = np.array((-a, -a, a / 3.0), dtype=float)
    point2 = make_point(array2) 
    return BRepPrimAPI_MakeBox(point1, point2)

def front_layer( side ):
    a = float(side / 2.0)
    array1 = -a * np.ones(3, dtype=float)
    point1 = make_point(array1)
    array2 = np.array((a, -a / 3.0, a), dtype=float)
    point2 = make_point(array2) 
    return BRepPrimAPI_MakeBox(point1, point2)

def back_layer( side ):
    a = float(side / 2.0)
    array1 = a * np.ones(3, dtype=float)
    point1 = make_point(array1)
    array2 = np.array((-a, a / 3.0, -a), dtype=float)
    point2 = make_point(array2) 
    return BRepPrimAPI_MakeBox(point1, point2)

def sphere_shape():
    diameter = 2.0
    origin = np.zeros(3, dtype=float)
    sphere = sphere_from_center_and_radius( origin,
                                            diameter / 2.0 )
    return sphere.Shape()
    
def bottom_cone():
    # create cone
    apex = np.array((-0.2, 0, 0), dtype=float)
    axis = np.array((1, 0, 0), dtype=float)
    height = 3.0
    radius1 = 0.0
    radius2 = 6.0

    cone = cone_from_point_height_directionvector_and_two_radii( apex,
                                                                 axis,
                                                                 height,
                                                                 radius1,
                                                                 radius2 )
    return cone.Shape()
   
def top_cone():
    # create cone
    apex = np.array((0.2, 0, 0), dtype=float)
    axis = np.array((-1, 0, 0), dtype=float)
    height = 3.0
    radius1 = 0.0
    radius2 = 6.0

    cone = cone_from_point_height_directionvector_and_two_radii( apex,
                                                                 axis,
                                                                 height,
                                                                 radius1,
                                                                 radius2 )
    return cone.Shape()
   
def front_cone():
    # create cone
    apex = np.array((0, 0.2, 0), dtype=float)
    axis = np.array((0, -1, 0), dtype=float)
    height = 3.0
    radius1 = 0.0
    radius2 = 6.0

    cone = cone_from_point_height_directionvector_and_two_radii( apex,
                                                                 axis,
                                                                 height,
                                                                 radius1,
                                                                 radius2 )
    return cone.Shape()

def back_cone():
    # create cone
    apex = np.array((0, -0.2, 0), dtype=float)
    axis = np.array((0, 1, 0), dtype=float)
    height = 3.0
    radius1 = 0.0
    radius2 = 6.0

    cone = cone_from_point_height_directionvector_and_two_radii( apex,
                                                                 axis,
                                                                 height,
                                                                 radius1,
                                                                 radius2 )
    return cone.Shape()

   
def left_cone():
    # create cone
    apex = np.array((0, 0, 0.2), dtype=float)
    axis = np.array((0, 0, -1), dtype=float)
    height = 3.0
    radius1 = 0.0
    radius2 = 6.0

    cone = cone_from_point_height_directionvector_and_two_radii( apex,
                                                                 axis,
                                                                 height,
                                                                 radius1,
                                                                 radius2 )
    return cone.Shape()

def right_cone():
    # create cone
    apex = np.array((0, 0, -0.2), dtype=float)
    axis = np.array((0, 0, 1), dtype=float)
    height = 3.0
    radius1 = 0.0
    radius2 = 6.0

    cone = cone_from_point_height_directionvector_and_two_radii( apex,
                                                                 axis,
                                                                 height,
                                                                 radius1,
                                                                 radius2 )
    return cone.Shape()

def cube_shape():
    side = 3.0
    return centered_cube_from_side(side).Shape()

def draw_sphere(event = None):
    display.DisplayColoredShape( sphere_shape() , 'RED' )

def draw_cone_1(event=None):
    display.DisplayColoredShape( bottom_cone() , 'GREEN' )

def draw_cone_2(event=None):
    display.DisplayColoredShape( top_cone() , 'GREEN' )

def draw_cube(event=None):
    display.DisplayColoredShape( cube_shape() , 'WHITE' )

def cut_bottom_cone():
    sphereShape = sphere_shape()
    coneShape = bottom_cone()
    cone_minus_sphere = BRepAlgoAPI_Cut(coneShape, sphereShape).Shape()
    return BRepAlgoAPI_Common(cone_minus_sphere, cube_shape()).Shape()

def cut_top_cone():
    sphereShape = sphere_shape()
    coneShape = top_cone()
    cone_minus_sphere = BRepAlgoAPI_Cut(coneShape, sphereShape).Shape()
    return BRepAlgoAPI_Common(cone_minus_sphere, cube_shape()).Shape()

def bottom_cone_layer():
    sphereShape = sphere_shape()
    coneShape = bottom_cone()
    layer = bottom_layer(3.0).Shape()
    step1 =  BRepAlgoAPI_Common(coneShape, cube_shape()).Shape()
    step2 =  BRepAlgoAPI_Fuse(step1, layer).Shape()
    return BRepAlgoAPI_Cut(step2, sphereShape).Shape()

def top_cone_layer():
    sphereShape = sphere_shape()
    coneShape = top_cone()
    layer = top_layer(3.0).Shape()
    step1 =  BRepAlgoAPI_Common(coneShape, cube_shape()).Shape()
    step2 =  BRepAlgoAPI_Fuse(step1, layer).Shape()
    return BRepAlgoAPI_Cut(step2, sphereShape).Shape()

def front_cone_layer():
    sphereShape = sphere_shape()
    coneShape = front_cone()
    layer = front_layer(3.0).Shape()
    step1 =  BRepAlgoAPI_Common(coneShape, cube_shape()).Shape()
    step2 =  BRepAlgoAPI_Fuse(step1, layer).Shape()
    return BRepAlgoAPI_Cut(step2, sphereShape).Shape()

def right_cone_layer():
    sphereShape = sphere_shape()
    coneShape = right_cone()
    layer = right_layer(3.0).Shape()
    step1 =  BRepAlgoAPI_Common(coneShape, cube_shape()).Shape()
    step2 =  BRepAlgoAPI_Fuse(step1, layer).Shape()
    return BRepAlgoAPI_Cut(step2, sphereShape).Shape()

def left_cone_layer():
    sphereShape = sphere_shape()
    coneShape = left_cone()
    layer = left_layer(3.0).Shape()
    step1 =  BRepAlgoAPI_Common(coneShape, cube_shape()).Shape()
    step2 =  BRepAlgoAPI_Fuse(step1, layer).Shape()
    return BRepAlgoAPI_Cut(step2, sphereShape).Shape()

def back_cone_layer():
    sphereShape = sphere_shape()
    coneShape = back_cone()
    layer = back_layer(3.0).Shape()
    step1 = BRepAlgoAPI_Common(coneShape, cube_shape()).Shape()
    step2 = BRepAlgoAPI_Fuse(step1, layer).Shape()
    return BRepAlgoAPI_Cut(step2, sphereShape).Shape()

def face_piece():
    back = back_cone_layer()
    left = left_cone_layer()
    right = right_cone_layer()
    top = top_cone_layer()
    bottom = bottom_cone_layer()
    step1 = BRepAlgoAPI_Cut(back, left).Shape()
    step2 = BRepAlgoAPI_Cut(step1, right).Shape()
    step3 = BRepAlgoAPI_Cut(step2, top).Shape()
    return BRepAlgoAPI_Cut(step3, bottom).Shape()

def edge_piece():
    back = back_cone_layer()
    left = left_cone_layer()
    top = top_cone_layer()
    bottom = bottom_cone_layer()
    step1 = BRepAlgoAPI_Common(back, left).Shape()
    step2 = BRepAlgoAPI_Cut(step1, top).Shape()
    return BRepAlgoAPI_Cut(step2, bottom).Shape()

def corner_piece():
    back = back_cone_layer()
    left = left_cone_layer()
    top = top_cone_layer()
    step = BRepAlgoAPI_Common(back, left).Shape()
    return BRepAlgoAPI_Common(step, top).Shape()

def draw_cut_bottom_cone(event=None):
    display.DisplayColoredShape( cut_bottom_cone() , 'BLACK' )

def draw_cut_top_cone(event=None):
    display.DisplayColoredShape( cut_top_cone() , 'BLACK' )

def draw_bottom_layer(event=None):
    display.DisplayColoredShape( bottom_cone_layer() , 'BLUE' )

def draw_top_layer(event=None):
    display.DisplayColoredShape( top_cone_layer() , 'BLUE' )

def draw_back_layer(event=None):
    display.DisplayColoredShape( front_cone_layer() , 'BLUE' )

def draw_front_layer(event=None):
    display.DisplayColoredShape( back_cone_layer() , 'BLUE' )

def draw_left_layer(event=None):
    display.DisplayColoredShape( left_cone_layer() , 'BLUE' )

def draw_right_layer(event=None):
    display.DisplayColoredShape( right_cone_layer() , 'BLUE' )

def draw_face_piece(event=None):
    display.DisplayColoredShape( face_piece() , 'CYAN' )

def draw_edge_piece(event=None):
    display.DisplayColoredShape( edge_piece() , 'GREEN' )

def draw_corner_piece(event=None):
    display.DisplayColoredShape( corner_piece() , 'ORANGE' )

def erase_all(event=None):
    display.EraseAll()

if __name__ == '__main__':
    toSTL(corner_piece(), 'corner.stl')
    toSTL(edge_piece(), 'edge.stl')
    toSTL(face_piece(), 'face.stl')

    display, start_display, add_menu, add_function_to_menu = \
        init_display()
    add_menu('Draw')
    add_function_to_menu('Draw', draw_sphere)
    add_function_to_menu('Draw', draw_cone_1)
    add_function_to_menu('Draw', draw_cone_2)
    add_function_to_menu('Draw', draw_cube)
    add_function_to_menu('Draw', draw_cut_bottom_cone)
    add_function_to_menu('Draw', draw_cut_top_cone)
    add_menu('Cubes')
    add_function_to_menu('Cubes', draw_face_piece)
    add_function_to_menu('Cubes', draw_edge_piece)
    add_function_to_menu('Cubes', draw_corner_piece)
    add_menu('Layers')
    add_function_to_menu('Layers', draw_bottom_layer)
    add_function_to_menu('Layers', draw_top_layer)
    add_function_to_menu('Layers', draw_front_layer)
    add_function_to_menu('Layers', draw_back_layer)
    add_function_to_menu('Layers', draw_left_layer)
    add_function_to_menu('Layers', draw_right_layer)
    add_menu('Erase')
    add_function_to_menu('Erase', erase_all)
    start_display()
