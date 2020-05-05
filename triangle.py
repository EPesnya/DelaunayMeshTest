import os
import sys
import numpy as np
import plotly.graph_objects as go

#############

n_p = 0
n_tri = 0
points = []
triangles = []
Delone = True
n_bad = 0

#############

f = open(os.path.join(sys.path[0], "strange_prism.vtk"),"r")
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

#################

def norm(v):
    return (v[0]**2 + v[1]**2 + v[2]**2)**0.5


def centr_circle(A, B, C):
    AC_centr = [0.5*(C[0] + A[0]), 0.5*(C[1] + A[1]), 0.5*(C[2] + A[2])]
    AB_centr = [0.5*(A[0] + B[0]), 0.5*(A[1] + B[1]), 0.5*(A[2] + B[2])]
    AB = [B[0] - A[0], B[1] - A[1], B[2] - A[2]]
    AC = [C[0] - A[0], C[1] - A[1], C[2] - A[2]]
    M = np.array([
        [AC[0], AC[1], AC[2]], [AB[0], AB[1], AB[2]],
        [AB[1]*AC[2] - AB[2]*AC[1], AB[2]*AC[0] - AB[0]*AC[2], AB[0]*AC[1] - AC[0]*AB[1]]
        ])
    v = np.array([
        sum(list(map(lambda a,b : a*b,AC,AC_centr))),
        sum(list(map(lambda a,b : a*b,AB,AB_centr))),
        A[0]*(AB[1]*AC[2] - AB[2]*AC[1]) + A[1]*(AB[2]*AC[0] - AB[0]*AC[2]) + A[2]*(AB[0]*AC[1] - AC[0]*AB[1])
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

    det_niz = np.array([
        [2*(S[0] - A[0]), 2*(S[1] - A[1]), 2*(S[2] - A[2])],
        [2*(S[0] - B[0]), 2*(S[1] - B[1]), 2*(S[2] - B[2])],
        [2*(S[0] - C[0]), 2*(S[1] - C[1]), 2*(S[2] - C[2])]
        ])

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

    for tri in triangles:
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
                zero_dist_index = np.where(abs_dist == min(abs_dist))
                zero_b_point = (b_points[zero_dist_index])[0]

                bad_dot.append(go.Scatter3d(
                    x = [zero_b_point[0]],
                    y = [zero_b_point[1]],
                    z = [zero_b_point[2]]
                    ))

                tri_plot.append(go.Scatter3d(
                    x = [a[0], b[0], c[0]],
                    y = [a[1], b[1], c[1]],
                    z = [a[2], b[2], c[2]],
                    mode='markers',
                    marker=dict(
                        color=[2, 2, 2],
                        colorscale='Viridis',
                        opacity=0.8
                        )
                    ))

                x_s = r*np.cos(u)*np.sin(v)
                y_s = r*np.sin(u)*np.sin(v)
                z_s = r*np.cos(v)
                sphere.append(go.Scatter3d(
                    x = x_s.flatten() + center_0[0],
                    y = y_s.flatten() + center_0[1],
                    z = z_s.flatten() + center_0[2],
                    mode='lines'
                    ))

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

            # Draw

            bad_dot.append(go.Scatter3d(
                x = [min_b_point[0]],
                y = [min_b_point[1]],
                z = [min_b_point[2]]
                ))

            tri_plot.append(go.Scatter3d(
                x = [a[0], b[0], c[0]],
                y = [a[1], b[1], c[1]],
                z = [a[2], b[2], c[2]],
                mode='markers',
                marker=dict(
                    color=[2, 2, 2],
                    colorscale='Viridis',
                    opacity=0.8
                    )
                ))

            x_s = r*np.cos(u)*np.sin(v)
            y_s = r*np.sin(u)*np.sin(v)
            z_s = r*np.cos(v)
            sphere.append(go.Scatter3d(
                x = x_s.flatten() + center_0[0],
                y = y_s.flatten() + center_0[1],
                z = z_s.flatten() + center_0[2],
                mode='lines'
                ))

            x_s = R*np.cos(u)*np.sin(v)
            y_s = R*np.sin(u)*np.sin(v)
            z_s = R*np.cos(v)
            sphere.append(go.Scatter3d(
                x = x_s.flatten() + center_1[0],
                y = y_s.flatten() + center_1[1],
                z = z_s.flatten() + center_1[2],
                mode='lines'
                ))

            for point in points:
                if norm(point - center_1) < R:
                    return False

    return True

############


points = np.array(points)
triangles = np.array(triangles)

print(check_mesh_delaunay(points, triangles))

tri_points = points[triangles]

# #extract the lists of x, y, z coordinates of the triangle vertices and connect them by a line
Xe = []
Ye = []
Ze = []
for T in tri_points:
    Xe.extend([T[k%3][0] for k in range(4)]+[ None])
    Ye.extend([T[k%3][1] for k in range(4)]+[ None])
    Ze.extend([T[k%3][2] for k in range(4)]+[ None])
       
# #define the trace for triangle sides
lines = go.Scatter3d(
                   x=Xe,
                   y=Ye,
                   z=Ze,
                   mode='lines'
                )


# surface = go.Mesh3d(
#     x=points[:,0],
#     y=points[:,1],
#     z=points[:,2],
#     i=triangles[:,0],
#     j=triangles[:,1],
#     k=triangles[:,2]
#     )

fig = go.Figure(data=[lines] + sphere + tri_plot + bad_dot)
fig.show()
