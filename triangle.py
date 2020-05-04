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

points = np.array(points)
triangles = np.array(triangles)

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
        [AB[1]*AC[2] - AB[2]*AC[1], AB[2]*AC[0] - AB[0]*AC[2], AB[0]*AC[1] - AC[0]*AB[1]]])
    v = np.array([
        sum(list(map(lambda a,b : a*b,AC,AC_centr))),
        sum(list(map(lambda a,b : a*b,AB,AB_centr))),
        A[0]*(AB[1]*AC[2] - AB[2]*AC[1]) + A[1]*(AB[2]*AC[0] - AB[0]*AC[2]) + A[2]*(AB[0]*AC[1] - AC[0]*AB[1])
        ])
    return np.linalg.solve(M, v)


def t_max(a, k, n):
    b = 0
    sum_of_sq_a_k = 0
    for i in range(len(a)):
        b += n[i]*(a[i] - k[i])
        n_in_sq = [i**2 for i in n]
        n_sum_sq = sum(n_in_sq)
    for i in range(len(a)):
        b += n[i]*(a[i] - k[i])
    for i in range(len(a)):
        sum_of_sq_a_k += (a[i] + k[i])**2
    D = b**2 - n_sum_sq * (-r_max**2 + sum_of_sq_a_k)
    print('D = {}'.format(D))
    t = (b + D**0.5) / n_sum_sq
    return t

################

b_min = np.array([min(points[:, 0]), min(points[:, 1]), min(points[:, 2])])
b_max = np.array([max(points[:, 0]), max(points[:, 1]), max(points[:, 2])])
r_max = norm(b_max - b_min)

u = np.linspace(0, 2*np.pi, 10)
v = np.linspace(0, np.pi, 10)
u, v = np.meshgrid(u, v)
sphere = []
one_tri = None

kostil = 0
################

for i in range(n_tri):
    tri = triangles[i]

    a = points[tri[0]]
    b = points[tri[1]]
    c = points[tri[2]]

    A = (b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1])
    B = (b[0] - a[0]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[0] - a[0])
    C = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])

    norma = norm([A, B, C])
    A = A / norma
    B = B / norma
    C = C / norma
    normal = np.array([A, B, C])

    center_0 = centr_circle(a, b, c)

    # ax.plot(
    #     [a[0], a[0] + normal[0]],
    #     [a[1], a[1] + normal[1]],
    #     [a[2], a[2] + normal[2]], 'r')

    # ax.plot(
    #     [a[0], b[0], c[0], a[0]],
    #     [a[1], b[1], c[1], a[1]],
    #     [a[2], b[2], c[2], a[2]], color='g')

    # if kostil == 0:
    #     one_tri = go.Scatter3d(
    #             x = [a[0], b[0], c[0]],
    #             y = [a[1], b[1], c[1]],
    #             z = [a[2], b[2], c[2]],
    #             mode='markers',
    #             marker=dict(
    #                 color=[1,2,3], # set color to an array/list of desired values
    #                 colorscale='Viridis', # choose a colorscale
    #                 opacity=0.8
    #                 )
    #         )


    good_triangle = False
    t_list = [0, 0.1, 1, 5, 10, 100]

    for t in t_list:
        center_1 = center_0 + t * normal
        # center_2 = center_0 - t * normal
        r = norm(a - center_1)
        good_sphere = True

        for j in range(n_p):
            if j not in tri:
                point = points[j]

                if norm(point - center_1) < r:
                    # print('tiangle: ' + str(i) + '; dot: ' + str(j))
                    # ax.plot(
                    #     [a[0], b[0], c[0], a[0]],
                    #     [a[1], b[1], c[1], a[1]],
                    #     [a[2], b[2], c[2], a[2]], color='g')
                    # ax.scatter(point[0], point[1], point[2], color='r')
                    # if kostil == 0:
                    #     x_s = r*np.cos(u)*np.sin(v)
                    #     y_s = r*np.sin(u)*np.sin(v)
                    #     z_s = r*np.cos(v)
                    #     sphere.append(go.Scatter3d(
                    #                 x = x_s.flatten() + center_1[0],
                    #                 y = y_s.flatten() + center_1[1],
                    #                 z = z_s.flatten() + center_1[2],
                    #                 mode='lines'
                    #             ))

                    good_sphere = False
                    break
            
        if good_sphere == True:
            good_triangle = True
            break

    if good_triangle == False:
        Delone = False
        n_bad += 1

        # # Display spheres
        r = norm(a - center_0)
        x_s = r*np.cos(u)*np.sin(v)
        y_s = r*np.sin(u)*np.sin(v)
        z_s = r*np.cos(v)

        # sphere.append(go.Scatter3d(
        #         x = x_s.flatten() + center_0[0],
        #         y = y_s.flatten() + center_0[1],
        #         z = z_s.flatten() + center_0[2],
        #         opacity=0.50
        #     ))

        kostil = 1

        
print(Delone)
print('n_bad = ' + str(n_bad) + '\n')

############


tri_points = points[triangles]

#extract the lists of x, y, z coordinates of the triangle vertices and connect them by a line
Xe = []
Ye = []
Ze = []
for T in tri_points:
    Xe.extend([T[k%3][0] for k in range(4)]+[ None])
    Ye.extend([T[k%3][1] for k in range(4)]+[ None])
    Ze.extend([T[k%3][2] for k in range(4)]+[ None])
       
#define the trace for triangle sides
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

fig = go.Figure(data=[lines] + sphere)


fig.show()

# plt.show()
