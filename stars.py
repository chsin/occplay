import operator
import cadmium
import itertools
import numpy as np

def tolerance():
    return 0.0001

def mag(n):
    return np.sqrt(sum(n **2))

def rot(solid, x = 0, y= 0, z = 0):
    # does two rotations:
    # 1. about z-axis
    # 2. about cross product of input vector and z-axis

    a = np.array([x, y, z]) / mag(np.array([x, y, z]))
    Zaxis = np.array([0.0, 0.0, 1.0], dtype = float)
    Caxis = np.cross(Zaxis, a)

    if (mag(Caxis) > tolerance()):
        c = Caxis / mag(Caxis)

        Zangle = 180.0 * np.arccos(c[0]) / np.pi
        solid.rotate(axis = list(Zaxis), angle = Zangle)

        Cangle = 180.0 * np.arccos(a[2]) / np.pi
        solid.rotate(axis = list(Caxis), angle = Cangle)

    return solid

def sum(seq):
    if len(seq) == 0:
        return 0 
    return reduce(operator.add, seq)

def intersect(seq):
    if len(seq) == 0:
        return 0
    return reduce(operator.mul, seq)

def pumpkin():
    spheres = [cadmium.Sphere(r = 1.1).translate(np.sin(2.0 * np.pi * r), 
                                                 np.cos(2.0 * np.pi * r), 
                                                 0) for r in np.linspace(0.0, 1.0, 12)]
    sum(spheres).toSTL("pumkin.stl")

def cubeVerts():
    # cube verts on a unit sphere
    points = [(1, 1, 1), (1, 1, -1), (1, -1, 1), (1, -1, -1), 
              (-1, 1, 1), (-1, 1, -1), (-1, -1, 1), (-1, -1, -1)]
    verts = [np.array(p) / mag(np.array(p)) for p in points]
    return verts

def tetrahedronVerts():
    # tetrahedron verts on a unit sphere centered at the origin
    psi = 1.0 / np.sqrt(2)
    points = [(1, 0, -psi), (-1, 0, -psi), (0, 1, psi), (0, -1, psi)]
    verts = [np.array(p) / mag(np.array(p)) for p in points]
    return verts

def icosahedronVerts():
    # icosahedron verts on a unit sphere centered at the origin
    phi = (1.0 + np.sqrt(5))/ 2.0
    points = [(0, 1, phi), (0, -1, phi), (0, 1, -phi), (0, -1, -phi),
              (1, phi, 0), (-1, phi, 0), (1, -phi, 0), (-1, -phi, 0),
              (phi, 0, 1), (phi, 0, -1), (-phi, 0, 1), (-phi, 0, -1)]
    verts = [np.array(p) / mag(np.array(p)) for p in points]
    return verts

def cube(side = 1.0):
    cube = cadmium.Box(x = side, y = side, z = side, center = True)
    cube.toSTL("cube.stl")
    return cube

def octahedron(side = 1.0):
    verts = cubeVerts()
    cubes = [rot(cadmium.Box(x = 1.5 * side, y = 1.5 * side, z = 1.5 * side, 
                             center = True), 
                 *tuple(v)).translate(*tuple(side * v)) 
             for v in verts]
    octa = (cadmium.Sphere(r = 0.7 * side, center = True) - sum(cubes))
    octa.toSTL("octa.stl")
    return octa

def tetrahedron(side = 1.0):
    verts = tetrahedronVerts()
    norms = closestNorm3(verts)

    cubes = [rot(cadmium.Box(x = 2.0 * side, y = 2.0 * side, z = 2.0 * side, 
                             center = True), 
                 *tuple(n)).translate(*tuple(n * (mag(n) + side) / mag(n))) 
             for n in norms]
    tetra = (cadmium.Sphere(r = side, center = True) - sum(cubes))
    tetra.toSTL("tetra.stl")
    return tetra

def distance(verts):
    d = []
    for i, j in itertools.product(range(len(verts)), range(len(verts))):
        if i < j:
            d.append(round(np.sqrt(sum((verts[i] - verts[j]) ** 2)), 4))
    return sorted( set(d) )

def NewVertsAtDistance(verts, d):
    newverts = []
    for i, j in itertools.product(range(len(verts)), range(len(verts))):
            if  (i < j and np.abs(mag(verts[i] - verts[j]) - d) < tolerance()):
                newverts.append((verts[i] + verts[j]) / mag(verts[i] + verts[j]))
    return newverts

def closestThree(verts, d1, d2 = 0):
    if not d2: 
        d2 = d1

    centers = []
    L = len(verts)
    for i, j, k in itertools.product(range(L), range(L), range(L)):
        if  ((i < j) and 
             ((np.abs(mag(verts[i] - verts[j]) - d1) < tolerance()) or 
              (np.abs(mag(verts[i] - verts[j]) - d2) < tolerance()))):
            if  ((j < k) and 
                 ((np.abs(mag(verts[i] - verts[k]) - d1) < tolerance()) or 
                  (np.abs(mag(verts[i] - verts[k]) - d2) < tolerance())) and 
                 ((np.abs(mag(verts[j] - verts[k]) - d1) < tolerance()) or 
                  (np.abs(mag(verts[j] - verts[k]) - d2) < tolerance()))):
                for a, b, c in ((i, j, k), (j, k, i), (k, i, j)):
                    center = verts[a] + verts[b] + 2.8 * verts[c]
                    center = center / np.sqrt(sum(center ** 2))
                    centers.append(tuple([round(center[0], 4), round(center[1], 4), round(center[2], 4)]))
    return centers


def closestNorm4(verts):
    normals = []
    d, d2 = distance(verts)[:2]
    for i, j, k in itertools.product(range(len(verts)), range(len(verts)), range(len(verts))):
        if ((i < j) and (i < k)):
            if (np.abs(mag(verts[i] - verts[j]) - d) < tolerance()):
                if ((np.abs(mag(verts[i] - verts[k]) - d) < tolerance())
                    and (np.abs(mag(verts[j] - verts[k]) - d2) < tolerance())):
                    norm = np.cross(verts[i] - verts[j], verts[i] - verts[k])
                    if (mag(norm) > tolerance()):
                        normals.append(norm / mag(norm))
    return normals


def closestNorm3(verts):
    d =  distance(verts)[0]
    normals = []
    for i, j, k in itertools.product(range(len(verts)), range(len(verts)), range(len(verts))):
        if ((i < j) and (j < k)):
            if (np.abs(mag(verts[i] - verts[j]) - d) < tolerance()):
                if ((np.abs(mag(verts[i] - verts[k]) - d) < tolerance())
                    and (np.abs(mag(verts[j] - verts[k]) - d) < tolerance())):
                    norm = (verts[i] + verts[j] + verts[k]) / 3.0
                    normals.append(norm)
    return normals


def golf():
    verts = icosahedronVerts()

    d = distance(verts)[0]
    newverts = NewVertsAtDistance(verts, d)
    verts += newverts

    d1, d2 = distance(verts)[:2]
    centers = closestThree(verts, d1, d2)

    for v in newverts:
        centers.append(tuple([round(v[0], 4), round(v[1], 4), round(v[2], 4)]))
   
    spheres = [cadmium.Sphere(r = 0.1, center = True).translate(*v) for v in centers]
    print "subtracting ", len(spheres), " dimples"
       
    (cadmium.Sphere(r = 1.0, center = True) - sum(spheres)).toSTL("golf.stl")

def icosahedron(side = 1.0):
    verts = icosahedronVerts()
    norms = closestNorm3(verts)

    cubes = [rot(cadmium.Box(x = side, 
                             y = side, 
                             z = side, 
                             center = True), 
                 *tuple(n)).translate(*tuple(n * side)) 
             for n in norms]
    sphere = cadmium.Sphere(r = side / 2.0, center = True)
    icosa = (sphere - sum(cubes))
    # icosa.toSTL("icosa.stl")
    return icosa

def icosastar(side = 1.0):
    verts = icosahedronVerts()

    cubes = [cadmium.Sphere(r = side, center = True).translate(*tuple(v * 1.7 * side)) for v in verts]

    icosa = sum(cubes)

    return icosa


def dodecahedron(side = 1.0):
    verts = icosahedronVerts()
    cubes = [rot(cadmium.Box(x = 1.5 * side, y = 1.5 * side, z = 1.5 * side, 
                             center = True), 
                 *tuple(v)).translate(*tuple(v * side)) 
             for v in verts]
    dodeca = (cadmium.Sphere(r = side, center = True) - sum(cubes))
    return dodeca

def dodecahedronMilledStar(side = 1.0):
    verts = icosahedronVerts()
    cylinders = []
    powers = np.arange(10)
    for p in powers:
        a = 0.13 * side * 0.74**p
        cylinders += [rot(cadmium.Cylinder(r = a, h = a, 
                                           center = True), 
                          *tuple(v)).translate(*tuple(v * a * 2)) 
                      for v in verts]
    
    dodeca = (dodecahedron(side) - sum(cylinders))
    return dodeca

if __name__=='__main__':
    s = 114.0
    # export to stl
    dodecahedronMilledStar(s).toSTL("dodeca.stl")
    # export to step
    dodecahedronMilledStar(s).toSTEP("dodeca.step")
