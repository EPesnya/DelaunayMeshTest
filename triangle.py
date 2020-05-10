import numpy as np
import iovtk

#############

epsilon_r = 1e-10
epsilon_zero_dist = 1e-2

points = []
triangles = []
points, triangles = iovtk.read_ug_tri('cone')

l_u, l_v = 30, 30
u = np.linspace(0, 2*np.pi, l_u)
v = np.linspace(0, np.pi, l_v)
u, v = np.meshgrid(u, v)

x_s = np.cos(u)*np.sin(v)
y_s = np.sin(u)*np.sin(v)
z_s = np.cos(v)

s_points = []
spheres_elem = []

#################

def norm(v):
    return np.linalg.norm(v)

def norm_2(v):
    return v[0]**2 + v[1]**2 + v[2]**2

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

def plane_from_points(a, b, c):
    A = (b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1])
    B = (b[2] - a[2]) * (c[0] - a[0]) - (b[0] - a[0]) * (c[2] - a[2])
    C = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
    D = - a[0] * A - a[1] * B - a[2] * C
    return A, B, C, D

def add_sphere(x, y, z):
    n_s_p = len(s_points)
    s_points.append([x[0, 0], y[0, 0], z[0, 0]])

    for i in range(1, l_v - 1):
        for j in range(l_u - 1):
            s_points.append([x[i, j], y[i, j], z[i, j]])

    s_points.append([x[-1, 0], y[-1, 0], z[-1, 0]])

    for i in range(l_u - 2):
        spheres_elem.append([n_s_p, n_s_p + i + 1, n_s_p + i + 2])
    spheres_elem.append([n_s_p, n_s_p + l_u - 1, n_s_p + 1])

    for i in range(l_v - 3):
        for j in range(l_u - 2):
            spheres_elem.append([
                n_s_p + 1 + i * (l_u - 1) + j,
                n_s_p + 1 + (i + 1) * (l_u - 1) + j,
                n_s_p + 1 + (i + 1) * (l_u - 1) + j + 1,
                n_s_p + 1 + i * (l_u - 1) + j + 1
                ])
        spheres_elem.append([
            n_s_p + 1 + i * (l_u - 1) + l_u - 2,
            n_s_p + 1 + (i + 1) * (l_u - 1) + l_u - 2,
            n_s_p + 1 + (i + 1) * (l_u - 1),
            n_s_p + 1 + i * (l_u - 1)
            ])

    for i in range(l_u - 2):
        spheres_elem.append([
            n_s_p + 1 + (l_v - 3) * (l_u - 1) + i,
            n_s_p + 1 + (l_v - 3) * (l_u - 1) + i + 1,
            n_s_p + 1 + (l_v - 2) * (l_u - 1)
            ])
    spheres_elem.append([
        n_s_p + 1 + (l_v - 3) * (l_u - 1) + l_u - 2,
        n_s_p + 1 + (l_v - 3) * (l_u - 1),
        n_s_p + 1 + (l_v - 2) * (l_u - 1)
        ])

def find_neighbors(triangles):
    n = len(triangles)
    neighbors = [[] for i in range(n)]

    for i in range(n - 1):
        edge_1 = [triangles[i, 0], triangles[i, 1]]
        edge_2 = [triangles[i, 0], triangles[i, 2]]
        edge_3 = [triangles[i, 1], triangles[i, 2]]

        for j in range(i + 1, n):
            if set(edge_1).issubset(set(triangles[j])) or set(edge_2).issubset(set(triangles[j])) or set(edge_3).issubset(set(triangles[j])):
                neighbors[i].append(j)
                neighbors[j].append(i)

    return np.array(neighbors)

################

def check_mesh_delaunay(points, triangles):

    print('Checking...')
    bad_triangles = np.zeros(len(triangles), dtype=np.int8)

    print('Finding neighbors...')
    neighbors = find_neighbors(triangles)
    print('Neighbors are found!')


    def rebuild(tri_n, bad_neighbor):
        tri = triangles[tri_n]
        common = list(set(tri) & set(triangles[bad_neighbor]))

        if len(common) == 2:
            t_1 = list(set(tri) - set(common)) + [common[0]] + list(set(triangles[bad_neighbor]) - set(common))
            t_2 = list(set(triangles[bad_neighbor]) - set(common)) + [common[1]] + list(set(tri) - set(common))

            bad_triangles[bad_neighbor] = 1
            bad_triangles[tri_n] = 1

            sides_1 = [tri_n] + [list(set(neighbors[bad_neighbor]) - set([tri_n]))[0]] + \
                [list(set(neighbors[tri_n]) - set([bad_neighbor]))[1]]
            sides_2 = [bad_neighbor] + [list(set(neighbors[tri_n]) - set([bad_neighbor]))[0]] + \
                [list(set(neighbors[bad_neighbor]) - set([tri_n]))[1]]

            for j in range(3):
                if neighbors[list(set(neighbors[bad_neighbor]) - set([tri_n]))[1]][j] == bad_neighbor:
                    neighbors[list(set(neighbors[bad_neighbor]) - set([tri_n]))[1]][j] = tri_n
                if neighbors[list(set(neighbors[tri_n]) - set([bad_neighbor]))[1]][j] == tri_n:
                    neighbors[list(set(neighbors[tri_n]) - set([bad_neighbor]))[1]][j] = bad_neighbor

            neighbors[tri_n] = np.array(sides_1)
            neighbors[bad_neighbor] = np.array(sides_2)
            triangles[bad_neighbor] = np.array(t_1)
            triangles[tri_n] = np.array(t_2)

    def check_triangle(i):
        tri = triangles[i]
        tri_n = i

        a = points[tri[0]]
        b = points[tri[1]]
        c = points[tri[2]]

        center_0 = centr_circle(a, b, c)
        r_2 = norm_2(a - center_0)
        r = r_2**0.5

        b_points = []

        for i in range(len(points)):
            if i not in tri:
                if norm_2(points[i] - center_0) - r_2 < -epsilon_r:
                    b_points.append(points[i])


        if len(b_points) > 0:
            b_points = np.array(b_points)
            bad_neighbor = None

            for b_point in b_points:
                for n in neighbors[tri_n]:
                    pnts = points[triangles[n]]
                    if np.array_equal(b_point, pnts[0]) or np.array_equal(b_point, pnts[1]) or np.array_equal(b_point, pnts[2]):
                        bad_neighbor = n
                        break

            A, B, C, D = plane_from_points(a, b, c)
            dist = (
                A * b_points[:, 0] + \
                B * b_points[:, 1] + \
                C * b_points[:, 2] + \
                D
                )
            abs_dist = abs(dist)

            if bad_neighbor != None:
                normal_t = np.array([A, B, C]) 
                a_1 = points[triangles[bad_neighbor][0]]
                b_1 = points[triangles[bad_neighbor][1]]
                c_1 = points[triangles[bad_neighbor][2]]
                A, B, C, D = plane_from_points(a_1, b_1, c_1)
                normal_n = np.array([A, B, C])
                cos_angle_2 = np.dot(normal_t, normal_n)**2 / norm_2(normal_n) / norm_2(normal_t)

                if cos_angle_2 > 3/4:
                    rebuild(tri_n, bad_neighbor)
                    return True

            if min(abs_dist) < epsilon_zero_dist:
                add_sphere(
                    r*x_s + center_0[0],
                    r*y_s + center_0[1],
                    r*z_s + center_0[2]
                    )
                return False

            # Проверяем что все точки лежат с одной стороны

            for i in range(1, len(b_points)):
                if np.sign(dist[i]) != np.sign(dist[0]):
                    add_sphere(
                        r*x_s + center_0[0],
                        r*y_s + center_0[1],
                        r*z_s + center_0[2]
                        )
                    return False

            # Ищем ближайшую

            dist_to_center = \
                (b_points[:, 0] - center_0[0])**2 + \
                (b_points[:, 1] - center_0[1])**2 + \
                (b_points[:, 2] - center_0[2])**2

            min_dist_index = np.where(dist_to_center == min(dist_to_center))
            min_b_point = (b_points[min_dist_index])[0]

            center_1 = centr_sphere(min_b_point, a, b, c)
            R_2 = norm_2(a - center_1)

            for point in points:
                if norm_2(point - center_1) - R_2 < -epsilon_r:
                    add_sphere(
                        r*x_s + center_0[0],
                        r*y_s + center_0[1],
                        r*z_s + center_0[2]
                        )
                    return False

        return True

    for i in range(len(triangles)):
        if bad_triangles[i] == 0:
            if check_triangle(i) == False:
                bad_triangles[i] = 2

    print('Checked!')
    return bad_triangles

#############

points = np.array(points)
triangles = np.array(triangles)

bad_triangles = check_mesh_delaunay(points, triangles)


iovtk.write_ug('output', points, triangles, bad_triangles)

if len(s_points) > 0:
    iovtk.write_ug('output_spheres', s_points, spheres_elem)