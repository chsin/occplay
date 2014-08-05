import sys
from OCC.Display.SimpleGui import *
from OCC.gp import *
from OCC.BRepPrimAPI import *
from OCC.BRepAlgoAPI import *
import numpy as np
from scipy import linalg

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

def corner_cube1( side ):
    a = float(side / 2.0);
    array1 = np.array((1, 1, 1), dtype=float)
    point1 = make_point(-a * array1)
    point2 = make_point((-a + side / 3.0) * array1)
    return BRepPrimAPI_MakeBox(point1, point2)

def corner_cube2( side ):
    a = float(side / 2.0);
    array2 = np.array((1, 1, -1), dtype=float)
    point1 = make_point(-a * array2)
    point2 = make_point((-a + side / 3.0) * array2)
    return BRepPrimAPI_MakeBox(point1, point2)

def corner_cube3( side ):
    a = float(side / 2.0);
    array3 = np.array((1, -1, 1), dtype=float)
    point1 = make_point(-a * array3)
    point2 = make_point((-a + side / 3.0) * array3)
    return BRepPrimAPI_MakeBox(point1, point2)

def corner_cube4( side ):
    a = float(side / 2.0);
    array4 = np.array((1, -1, -1), dtype=float)
    point1 = make_point(-a * array4)
    point2 = make_point((-a + side / 3.0) * array4)
    return BRepPrimAPI_MakeBox(point1, point2)

def sphere_shape():
    diameter = 1.5
    origin = np.zeros(3, dtype=float)
    sphere = sphere_from_center_and_radius( origin,
                                            diameter / 2.0 )
    return sphere.Shape()
    
def cone_shape():
    # create cone
    apex = np.zeros(3, dtype=float)
    axis = np.array((1, 0, 0), dtype=float)
    height = 2.0
    radius1 = 0.0
    radius2 = 4.0

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

def draw_cone(event=None):
    display.DisplayColoredShape( cone_shape() , 'GREEN' )

def draw_cube(event=None):
    display.DisplayColoredShape( cube_shape() , 'WHITE' )

def cut_cone_shape():
    sphereShape = sphere_shape()
    coneShape = cone_shape()
    cone_minus_sphere = BRepAlgoAPI_Cut(coneShape, sphereShape).Shape()
    
    return BRepAlgoAPI_Common(cone_minus_sphere, cube_shape()).Shape()

def corner_piece_shape1():
    cube = corner_cube1( 3.0 ).Shape()
    return BRepAlgoAPI_Cut(cube, sphere_shape()).Shape()

def corner_piece_shape2():
    cube = corner_cube2( 3.0 ).Shape()
    return BRepAlgoAPI_Cut(cube, sphere_shape()).Shape()

def corner_piece_shape3():
    cube = corner_cube3( 3.0 ).Shape()
    return BRepAlgoAPI_Cut(cube, sphere_shape()).Shape()

def corner_piece_shape4():
    cube = corner_cube4( 3.0 ).Shape()
    return BRepAlgoAPI_Cut(cube, sphere_shape()).Shape()

def draw_corner_piece_1():
    display.DisplayColoredShape( corner_piece_shape1() , 'BLACK' )

def draw_corner_piece_2():
    display.DisplayColoredShape( corner_piece_shape2() , 'BLACK' )

def draw_corner_piece_3():
    display.DisplayColoredShape( corner_piece_shape3() , 'BLACK' )

def draw_corner_piece_4():
    display.DisplayColoredShape( corner_piece_shape4() , 'BLACK' )

def draw_cut_cone(event=None):
    display.DisplayColoredShape( cut_cone_shape() , 'BLACK' )

def erase_all(event=None):
    display.EraseAll()

if __name__ == '__main__':
    display, start_display, add_menu, add_function_to_menu = \
        init_display()
    add_menu('Draw')
    add_function_to_menu('Draw', draw_sphere)
    add_function_to_menu('Draw', draw_cone)
    add_function_to_menu('Draw', draw_cube)
    add_function_to_menu('Draw', draw_cut_cone)
    add_menu('Cubes')
    add_function_to_menu('Cubes', draw_corner_piece_1)
    add_function_to_menu('Cubes', draw_corner_piece_2)
    add_function_to_menu('Cubes', draw_corner_piece_3)
    add_function_to_menu('Cubes', draw_corner_piece_4)
    add_menu('Erase')
    add_function_to_menu('Erase', erase_all)
    start_display()
