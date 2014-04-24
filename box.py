import operator
import cadmium
import itertools
import numpy as np

def rot(solid, x = 0, y= 0, z = 0):
    # holly fuck angle in Degrees?! 
    # who in the world does geometry and uses degrees
    # done bitching
    a = np.array([x, y, z])
    angle = 180 * a / np.sqrt(sum(np.array(a) ** 2)) / np.pi
    axis = np.array([cadmium.X_axis, cadmium.Y_axis, cadmium.Z_axis])
    for i in range(3):      
        solid.rotate(axis = axis[i], angle = angle[i])
    return solid

def sum(seq):
    if len(seq) == 0:
        return 0 
    return reduce(operator.add, seq)

def pumpkin():
    spheres = [cadmium.Sphere(r = 1.1).translate(np.sin(2.0 * np.pi * r), 
                                                 np.cos(2.0 * np.pi * r), 
                                                 0) for r in np.linspace(0.0, 1.0, 12)]
    sum(spheres).toSTL("pumkin.stl")

def cubeVerts():
    # cube verts on a unit sphere
    points = [(1, 1, 1), (1, 1, -1), (1, -1, 1), (1, -1, -1), 
              (-1, 1, 1), (-1, 1, -1), (-1, -1, 1), (-1, -1, -1)]
    verts = [np.array(p) / np.sqrt(sum(np.array(p) ** 2)) for p in points]
    return verts

def cube():
    verts = cubeVerts()
    side = 1.0
    cubes = [rot(cadmium.Box(x = side, y = side, z = side, 
                             center = True).translate(*v),
                 *v) for v in verts]

    sum(cubes).toSTL("cubes.stl")
    (cadmium.Sphere(r = 0.5, center = True) - sum(cubes)).toSTL("icosa.stl")

def icosahedronVerts():
    # icosahedron verts on a unit sphere
    phi = (1.0 + np.sqrt(5))/ 2.0
    points = [(0, 1, phi), (0, -1, phi), (0, 1, -phi), (0, -1, -phi),
              (1, phi, 0), (-1, phi, 0), (1, -phi, 0), (-1, -phi, 0),
              (phi, 0, 1), (phi, 0, -1), (-phi, 0, 1), (-phi, 0, -1)]
    verts = [np.array(p) / np.sqrt(sum(np.array(p) ** 2)) for p in points]
    return verts

def icosahedron():
    verts = icosahedronVerts()
    side = 1.0
    # cubes = [rot(cadmium.Box(x = a, y = a, z = a, 
    #                          center = True).translate(*v), 
    #              *v) for v in verts]
    cubes = []
    for v in verts:
        a = np.array(np.arccos(v), dtype = np.float)
        angle = 180 * a / np.sqrt(sum(np.array(a) ** 2)) / np.pi
        
        bok = cadmium.Box(x = side, y = side, z = side, center = True)
        bok.translate(*v)
        bok.rotate(axis = cadmium.X_axis, angle = angle[0])
        bok.rotate(axis = cadmium.Y_axis, angle = angle[0])
        bok.rotate(axis = cadmium.Z_axis, angle = angle[1])
        cubes.append(bok)
        

    (cadmium.Sphere(r = 1.0, center = True) - sum(cubes)).toSTL("icosa.stl")

def distance(verts):
    d = []
    for i, j in itertools.product(range(len(verts)), range(len(verts))):
        if i < j:
            d.append(round(sum((verts[i] - verts[j]) ** 2), 4))
    return sorted( set(d) )

def NewVertsAtDistance(verts, d):
    newverts = []
    for i, j in itertools.product(range(len(verts)), range(len(verts))):
            if  (i < j and np.abs(sum((verts[i] - verts[j]) ** 2) - d) < tolerance):
                newverts.append((verts[i] + verts[j]) / np.sqrt(sum((verts[i] + verts[j]) ** 2)))
    return newverts

def golf():
    verts = icosahedronVerts()

    d = distance(verts)[0]
    newverts = NewVertsAtDistance(verts, d)
    verts += newverts

    d1, d2 = distance(verts)[:2]
    centers = []
    for i, j, k in itertools.product(range(len(verts)), range(len(verts)), range(len(verts))):
        if  ((i < j) and ((np.abs(sum((verts[i] - verts[j]) ** 2) - d1) < tolerance) or (np.abs(sum((verts[i] - verts[j]) ** 2) - d2) < tolerance))):
            if  ((j < k) and ((np.abs(sum((verts[i] - verts[k]) ** 2) - d1) < tolerance) or (np.abs(sum((verts[i] - verts[k]) ** 2) - d2) < tolerance))
                  and ((np.abs(sum((verts[j] - verts[k]) ** 2) - d1) < tolerance) or (np.abs(sum((verts[j] - verts[k]) ** 2) - d2) < tolerance))):
                for a, b, c in ((i, j, k), (j, k, i), (k, i, j)):
                    center = verts[a] + verts[b] + 2.8 * verts[c]
                    center = center / np.sqrt(sum(center ** 2))
                    centers.append(tuple([round(center[0], 4), round(center[1], 4), round(center[2], 4)]))

    for v in newverts:
        centers.append(tuple([round(v[0], 4), round(v[1], 4), round(v[2], 4)]))
   
    spheres = [cadmium.Sphere(r = 0.1, center = True).translate(*v) for v in centers]
    print "subtracting ", len(spheres), " dimples"
       
    (cadmium.Sphere(r = 1.0, center = True) - sum(spheres)).toSTL("golf.stl")

if __name__=='__main__':
    tolerance = 0.0001
    # pumpkin()
    # golf()
    cube()
