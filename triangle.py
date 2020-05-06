import os
import sys
import numpy as np

#############

n_p = 0
n_tri = 0
points = []
triangles = []

#############

print('Reading file...')
f = open(os.path.join(sys.path[0], "cone.vtk"),"r")
fl = f.readlines()

p_line = 0
t_line = 0
for i in range(len(fl)):
    x = fl[i]
    if 'POINTS' in x:
        n_p = int(x.split(' ')[1])
        p_line = i + 1
    if 'CELLS' in x:
        n_tri = int(x.split(' ')[1])
        t_line = i + 1

for i in range(p_line, p_line + n_p):
    x = fl[i]
    ps = x.split(' ')
    points.append([float(ps[0]), float(ps[1]), float(ps[2])])

for i in range(t_line, t_line + n_tri):
    x = fl[i]
    ps = x.split(' ')
    triangles.append([int(ps[1]), int(ps[2]), int(ps[3])])

f.close()
print('Reading done!')

#################

def norm(v):
    return np.linalg.norm(v)

def centr_circle(A, B, C):
    AC_center = 0.5 * (A + C)
    AB_center = 0.5 * (A + B)
    AB = B - A
    AC = C - A
    M = np.array([AC, AB, np.cross(AB, AC)])
    v = np.array([
        sum(AC * AC_center),
        sum(AB * AB_center),
        sum(A * M[2])
        ])
    return np.linalg.solve(M, v)

def centr_sphere(A, B, C, S):

    x0_ver = np.array([
        [S[0]**2 - A[0]**2 + S[1]**2 - A[1]**2 + S[2]**2 - A[2]**2, 2*(S[1] - A[1]), 2*(S[2] - A[2])],
        [S[0]**2 - B[0]**2 + S[1]**2 - B[1]**2 + S[2]**2 - B[2]**2, 2*(S[1] - B[1]), 2*(S[2] - B[2])],
        [S[0]**2 - C[0]**2 + S[1]**2 - C[1]**2 + S[2]**2 - C[2]**2, 2*(S[1] - C[1]), 2*(S[2] - C[2])]
        ])

    y0_ver = np.array([
        [2*(S[0] - A[0]), S[0]**2 - A[0]**2 + S[1]**2 - A[1]**2 + S[2]**2 - A[2]**2, 2*(S[2] - A[2])],
        [2*(S[0] - B[0]), S[0]**2 - B[0]**2 + S[1]**2 - B[1]**2 + S[2]**2 - B[2]**2, 2*(S[2] - B[2])],
        [2*(S[0] - C[0]), S[0]**2 - C[0]**2 + S[1]**2 - C[1]**2 + S[2]**2 - C[2]**2, 2*(S[2] - C[2])]
        ])

    z0_ver = np.array([
        [2*(S[0] - A[0]), 2*(S[1] - A[1]), S[0]**2 - A[0]**2 + S[1]**2 - A[1]**2 + S[2]**2 - A[2]**2],
        [2*(S[0] - B[0]), 2*(S[1] - B[1]), S[0]**2 - B[0]**2 + S[1]**2 - B[1]**2 + S[2]**2 - B[2]**2],
        [2*(S[0] - C[0]), 2*(S[1] - C[1]), S[0]**2 - C[0]**2 + S[1]**2 - C[1]**2 + S[2]**2 - C[2]**2]
        ])

    det_niz = np.array([2*(S - A), 2*(S - B), 2*(S - C)])

    d = np.linalg.det(det_niz)
    x = np.linalg.det(x0_ver) / d
    y = np.linalg.det(y0_ver) / d
    z = np.linalg.det(z0_ver) / d
    return x, y, z 

################

u = np.linspace(0, 2*np.pi, 10)
v = np.linspace(0, np.pi, 10)
u, v = np.meshgrid(u, v)
sphere = []
tri_plot = []
bad_dot = []

################

def check_mesh_delaunay(points, triangles):

    bad_triangles = np.zeros(len(triangles), dtype=np.int8)

    def check_triangle(i):
        tri = triangles[i]

        a = points[tri[0]]
        b = points[tri[1]]
        c = points[tri[2]]

        center_0 = centr_circle(a, b, c)
        r = norm(a - center_0)

        b_points = []

        for i in range(len(points)):
            if i not in tri:
                if norm(points[i] - center_0) - r < -0.01:
                    b_points.append(points[i])

        if len(b_points) > 0:
            b_points = np.array(b_points)

            A = (b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1])
            B = (b[2] - a[2]) * (c[0] - a[0]) - (b[0] - a[0]) * (c[2] - a[2])
            C = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
            D = - a[0] * A - a[1] * B - a[2] * C

            norm_ABC = norm([A, B, C])

            dist = (
                A * b_points[:, 0] + \
                B * b_points[:, 1] + \
                C * b_points[:, 2] + D
                ) / norm_ABC

            # Проверяем не лежат ли 4 точки в одной плоскости

            abs_dist = abs(dist)

            if min(abs_dist) < 0.01:
                return False

            # Проверяем что все точки лежат с одной стороны

            for i in range(1, len(b_points)):
                if np.sign(dist[i]) != np.sign(dist[0]):
                    return False

            # Ищем ближайшую

            dist_to_center = \
                (b_points[:, 0] - center_0[0])**2 + \
                (b_points[:, 1] - center_0[1])**2 + \
                (b_points[:, 2] - center_0[2])**2

            min_dist_index = np.where(dist_to_center == min(dist_to_center))
            min_b_point = (b_points[min_dist_index])[0]

            center_1 = centr_sphere(min_b_point, a, b, c)
            R = norm(a - center_1)

            for point in points:
                if norm(point - center_1) - R < -0.01:
                    return False

        return True

    for i in range(len(triangles)):
        if check_triangle(i) == False:
            bad_triangles[i] = 1

    return bad_triangles

#############

points = np.array(points)
triangles = np.array(triangles)
bad_triangles = check_mesh_delaunay(points, triangles)

#############

print('Writing file...')
f = open(os.path.join(sys.path[0], "output.vtk"),"w")

f.write(
"""# vtk DataFile Version 2.0
Delaunay check output
ASCII
DATASET UNSTRUCTURED_GRID
"""
)

f.write('POINTS {} float\n'.format(n_p))
for i in range(n_p):
    f.write('{} {} {}\n'.format(points[i, 0], points[i, 1], points[i, 2]))

f.write('CELLS {} {}\n'.format(n_tri, 4 * n_tri))
for i in range(n_tri):
    f.write('3 {} {} {}\n'.format(triangles[i, 0], triangles[i, 1], triangles[i, 2]))

f.write('CELL_TYPES {}\n'.format(n_tri))
for i in range(n_tri):
    f.write('5\n')

f.write(
"""CELL_DATA {}
SCALARS scalars int
LOOKUP_TABLE default
""".format(n_tri)
)
for i in range(n_tri):
    f.write(str(bad_triangles[i]) + '\n')

f.close()
print('Writing done!')

###############
