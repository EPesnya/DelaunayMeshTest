import os
import sys

def read_ug_tri(filename):
    points = []
    triangles = []
    n_p = 0
    n_tri = 0
    filename += '.vtk'

    print('Reading file ' + filename)
    f = open(os.path.join(sys.path[0], filename),"r")
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

    return points, triangles


def write_ug(filename, points, elements, data = []):
    n_p = len(points)
    n_e = len(elements)
    element_length = [len(elements[i]) for i in range(n_e)]
    filename += '.vtk'

    print('Writing model ' + filename)
    f = open(os.path.join(sys.path[0], filename),"w")

    f.write(
        "# vtk DataFile Version 2.0\n" +
        "Delaunay check output\n" +
        "ASCII\n" +
        "DATASET UNSTRUCTURED_GRID\n"
        )

    f.write('POINTS {} float\n'.format(n_p))
    for i in range(n_p):
        f.write('{:.6f} {:.6f} {:.6f}\n'.format(
            points[i][0],
            points[i][1],
            points[i][2]))

    f.write('CELLS {} {}\n'.format(n_e, n_e + sum(element_length)))
    for i in range(n_e):
        f.write('{} {}\n'.format(
            element_length[i],
            ' '.join(map(str, elements[i]))
            ))

    f.write('CELL_TYPES {}\n'.format(n_e))
    for i in range(n_e):
        if element_length[i] == 3:
            f.write('5\n')
        else:
            f.write('9\n')

    if len(data) > 0:
        f.write(
        """CELL_DATA {}
        SCALARS scalars int
        LOOKUP_TABLE default
        """.format(n_e)
        )
        for i in range(n_e):
            f.write(str(data[i]) + '\n')

    f.close()
    print('Writing done!')