from packs import defpaths
import numpy as np
from packs.examples.benchmarks_biphasic.layers.generate_coarse_points import MoabMesh, Point, Line, Polygon
from typing import Sequence
import os
from pymoab import core, types, rng, topo_util



def write_points(points: Sequence[Point], lengths: list, points_ids_of_lenghts):

    text = []

    for point in points:
        id_lenght = None
        
        for i, list_points in enumerate(points_ids_of_lenghts):
            if point.id in list_points:
                id_lenght = i
                break

        idp = point.id + 1
        coord = point.coord
        coord_text = [str(value) for value in coord]
        coord_text += [str(lengths[id_lenght])]
        coord_text = ', '.join(coord_text)
        coord_text = '{' + coord_text + '};'
        line_text = 'Point(' + str(idp) + ')=' + coord_text + '\n'
        text.append(line_text)

        # content = 'Point(' + str(idp) + ')={' + 
    return text

def write_lines(lines: np.ndarray):

    n_lines = 1

    text = []
    for line in lines:
        line_text = 'Line(' + str(n_lines) + ')={' + str(line[0]) + ',' + str(line[1]) + '};' + '\n'
        text.append(line_text)
        n_lines += 1

    return text
    

def write_file_geo(filename, text_points, text_lines, text_curves):

    with open(filename, 'w') as f:
        for text in text_points:
            f.write(text)
        
        f.write('\n')
        
        for text in text_lines:
            f.write(text)
        
        f.write('\n')

        for text in text_curves:
            f.write(text)




def write_curve_loop(curves_loop):
    text = []
    n_curves = 1
    for curve in curves_loop:
        curve_str = [str(value) for value in curve]
        curve_str = ', '.join(curve_str)
        curve_str = '{' + curve_str + '};' 
        line = 'Curve Loop(' + str(n_curves) + ')=' + curve_str + '\n'
        text.append(line)
        line2 = 'Plane Surface(' + str(n_curves) + ')={' + str(n_curves) + '};' + '\n'
        text.append(line2)
        n_curves += 1
    
    return text

def write_coarse_msh_file(filename, elements, nodes_centroids):
    text = ['$MeshFormat\n']
    text.append('2.0 0 8\n')
    text.append('$EndMeshFormat\n')
    text.append('$Nodes\n')

    n_nodes = len(nodes_centroids)

    text.append(str(n_nodes) + '\n')

    for i, cnode in enumerate(nodes_centroids):
        coords_str = [str(val) for val in cnode]
        coords_str = ' '.join(coords_str)
        coords_str += '\n'
        linepoint = str(i+1) + ' ' + coords_str
        text.append(linepoint)
    
    text.append('$EndNodes\n')
    text.append('$Elements\n')

    n_elements = len(elements)

    text.append(str(n_elements) + '\n')

    n_points_elements = np.zeros(elements.shape[0], dtype=int)

    for i, el in enumerate(elements):
        n_points_elements[i] = len(el)
    
    test3 = n_points_elements == 3
    test4 = n_points_elements == 4

    triangles = elements[test3]
    quadrangles = elements[test4]

    type_triangle = '2'
    type_quadrangle = '3'

    n_el = 1

    for tri in triangles:
        # points_ids = tri - 1
        # coords = nodes_centroids[points_ids]
        coords_str = [str(val) for val in tri]
        coords_str = ' '.join(coords_str)
        coords_str += '\n'
        elem_str = str(n_el) + ' 2 0 ' + coords_str
        text.append(elem_str)
        n_el += 1

    for tri in quadrangles:
        # points_ids = tri - 1
        # coords = nodes_centroids[points_ids]
        coords_str = [str(val) for val in tri]
        coords_str = ' '.join(coords_str)
        coords_str += '\n'
        elem_str = str(n_el) + ' 3 0 ' + coords_str
        text.append(elem_str)
        n_el += 1
    
    text.append('$EndElements\n')

    with open(filename, 'w') as f:
        for line in text:
            f.write(line)
    
def write_coarse_mesh_pymoab(elements, points_list, nodes_centroids):

    # points_list = np.array([Point(coord[0], coord[1], coord[2]) for coord in nodes_centroids])
    polygons: Sequence[Polygon] = np.array([Polygon(points_list[np.array(indexes)-1]) for indexes in elements])

    for i, polygon in enumerate(polygons):
        polygon.sort_by_xy_plane()

    polygons_n_points = np.array([poly.npoints for poly in polygons])

    polygons3 = polygons[polygons_n_points==3]
    polygons4 = polygons[polygons_n_points==4]

    moab = MoabMesh()
    moab.create_vertexes(nodes_centroids)

    quads = [moab.verts[polygon.points_ids] for polygon in polygons4]
    tris = [moab.verts[polygon.points_ids] for polygon in polygons3]

    moab.create_tri_elements(tris)

    for i, quad in enumerate(quads):
        moab.create_quad_element(quad)
        
    moab.export_mesh('het_coarse2.vtk')
    moab.export_mesh('het_coarse2.msh')

    text = []
    n_el = 12

    for i, quad in enumerate(quads):
        
        poly: Polygon = polygons4[i]
        points_ids = poly.points_ids + 1

        coords_str = [str(val) for val in points_ids]
        coords_str = ' '.join(coords_str)
        coords_str += '\n'
        elem_str = str(n_el) + ' 3 0 ' + coords_str
        text.append(elem_str)
        n_el -= 1
        if n_el == 0:
            n_el -= 1
    
    file_text = os.path.join(defpaths.mesh, defpaths.biphasic_folder, 'correct_elements_het.txt')
    with open(file_text, 'w') as f:
        for line in text:
            f.write(line)
        





           





def run():

    filename = os.path.join(defpaths.mesh ,defpaths.biphasic_folder, 'coarse1_het.geo')
    filename_coarse_msh = os.path.join(defpaths.mesh ,defpaths.biphasic_folder, 'coarse1_het.msh')


    Lx = 1.0
    Ly = 1.0


    ## reta 1
    r1 = Line(0, 0.5, 0.5, 0)

    k1 = 1/6
    k2 = 1/3

    p1 = [0, 0]
    p2 = [k1, 0]
    p3 = [k1+k2, 0]
    p4 = [k1+2*k2, 0]
    p5 = [k1+2*k2 + k1, 0]
    p6 = [0, k1]
    p7 = [k2, k1]
    p8 = [2*k2, k1]
    p9 = [3*k2, k1]
    p10 = [p8[0] + Lx/6, 0]
    p11 = [p9[0] + Lx/6, p9[1]]
    p12 = [p10[0] + Lx/6, 0]
    p13 = [p11[0] + Lx/6, p11[1]]
    p14 = [p1[0], Ly/6]
    p15 = [p14[0] + Lx/6, p14[1]]
    p16 = [p15[0] + Lx/6, p15[1]]
    p17 = [p16[0] + Lx/6, p16[1]]
    p18 = [p17[0] + Lx/6, p16[1]]
    p19 = [p18[0] + Lx/6, p16[1]]


    points = np.array([p1, p2, p3, p4, p5, p6, p7, p8, p9])

    newpoints1 = points[[1, 2, 3]].copy()
    newpoints1[:, 1] = newpoints1[:, 1] + k2
    points = np.vstack([points, newpoints1])

    newpoints2 = points[[5, 6, 7, 8]].copy()
    newpoints2[:, 1] = newpoints2[:, 1] + k2
    points = np.vstack([points, newpoints2])

    newpoints3 = newpoints1.copy()
    newpoints3[:, 1] = newpoints3[:, 1] + k2
    points = np.vstack([points, newpoints3])

    newpoints4 = newpoints2.copy()
    newpoints4[:, 1] = newpoints4[:, 1] + k2
    points = np.vstack([points, newpoints4])

    newpoints5 = points[[0, 1, 2, 3, 4]].copy()
    newpoints5[:, 1] = newpoints5[:, 1] + 1
    points = np.vstack([points, newpoints5])

    zeros_append = np.zeros((points.shape[0], 1))
    points = np.hstack([points, zeros_append])


    lengths = [0.015]

    points_list = np.array([Point(coord[0], coord[1], coord[2]) for coord in points])
    all_points_ids = np.array([point.id for point in points_list])
    points_ids_lenght1 = np.arange(points.shape[0])

    points_ids_of_lenghts = np.array([
        points_ids_lenght1
    ], dtype='O')



    # lines = np.array([[points_list[0].id, points_list[1].id], [points_list[1].id, points_list[2].id], [points_list[2].id, points_list[0].id]])

    lines = np.array([
        [1, 2],
        [2, 6],
        [6, 1],
        [6, 10],
        [10, 13],
        [13, 6],
        [13, 17],
        [17, 20],
        [20, 13],
        [20, 25],
        [25, 24],
        [24, 20],
        [25, 21],
        [21, 26],
        [26, 25],
        [26, 22],
        [22, 27],
        [27, 26],
        [27, 23],
        [23, 28],
        [28, 27],
        [23, 19],
        [19, 16],
        [16, 23],
        [16, 12],
        [12, 9],
        [9, 16],
        [4, 5],
        [5, 9],
        [9, 4],
        [3, 4],
        [4, 8],
        [8, 3],
        [2, 3],
        [3, 7],
        [7, 2],
        [7, 10],
        [8, 11],
        [11, 7],
        [12, 8],
        [19, 15],
        [15, 12],
        [21, 18],
        [18, 22],
        [17, 21],
        [22, 19],
        [18, 15],
        [17, 14],
        [14, 18],
        [10, 14],
        [11, 14],
        [11, 15]

    ])
    
    curve_loop = np.array([ 
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9],
        [10, 11, 12],
        [13, 14, 15],
        [16, 17, 18],
        [19, 20, 21],
        [22, 23, 24],
        [25, 26, 27],
        [28, 29, 30],
        [31, 32, 33],
        [34, 35, 36],
        [37, -4, -2, -36],
        [38, 39, -35, -33],
        [40, -32, -30, -26],
        [41, 42, -25, -23],
        [43, 44, -16, -14],
        [45, -13, -10, -8],
        [46, -22, -19, -17],
        [47, -41, -46, -44],
        [48, 49, -43, -45],
        [50, -48, -7, -5],
        [51, -50, -37, -39],
        [52, -47, -49, -51],
        [-42, -52, -38, -40]
    ], dtype='O')

    elements2 = []
    for curve in curve_loop:
        cv = np.array(curve)
        cv = np.absolute(cv) - 1
        lines_cv = lines[cv]
        points_lines = np.unique(np.concatenate(lines_cv))
        elements2.append(points_lines)
    
    elements2 = np.array(elements2, dtype='O')

    # elements = np.array([
        
    # ], dtype='O')

    text_points = write_points(points_list, lengths, points_ids_of_lenghts)
    text_lines = write_lines(lines)
    text_curves = write_curve_loop(curve_loop)

    # write_file_geo(filename, text_points, text_lines, text_curves)
    # write_coarse_msh_file(filename_coarse_msh, elements, points)

    write_coarse_mesh_pymoab(elements2, points_list, points)

