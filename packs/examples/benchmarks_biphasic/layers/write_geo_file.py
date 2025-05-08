from packs import defpaths
import numpy as np
from packs.examples.benchmarks_biphasic.layers.generate_coarse_points import Point, Line
from typing import Sequence
import os


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

def run():

    filename = os.path.join(defpaths.mesh ,defpaths.layers_folder, 'coarse3.geo')

    Lx = 1.0
    Ly = 1.0


    ## reta 1
    r1 = Line(0, 0.5, 0.5, 0)

    p1 = [0, 0]
    p2 = [Lx/6, 0]
    p3 = [Lx/12, Ly/12]
    p4 = [2*Lx/6, 0]
    p5 = [p3[0] + Lx/6, Ly/12]
    p6 = [3*Lx/6, 0]
    p7 = [p5[0]+Lx/6, p5[1]]
    p8 = [p6[0] + Lx/6, 0]
    p9 = [p7[0] + Lx/6, p7[1]]
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


    points = np.array([p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14,
                       p15, p16, p17, p18, p19])
    
    newpoints1  = points[13:18].copy()
    newpoints1[:, 1] = newpoints1[:, 1] + Ly/12
    newpoints1[:, 0] = newpoints1[:, 0] + Lx/12

    points = np.vstack([points, newpoints1])

    newpoints2 = newpoints1.copy()
    newpoints2[:, 0] = newpoints2[:, 0] - Lx/12
    newpoints2[:, 1] = newpoints2[:, 1] + Ly/12
    points = np.vstack([points, newpoints2])

    newpoints3 = newpoints1[0:4]
    newpoints3[:, 1] = newpoints3[:, 1] + Ly/6
    points = np.vstack([points, newpoints3])

    points = np.hstack([points, np.zeros((points.shape[0], 1))])

    lengths = [0.02]

    points_list = np.array([Point(coord[0], coord[1], coord[2]) for coord in points])
    points_ids_of_lenghts = np.array([
        [point.id for point in points_list]
    ])

    # lines = np.array([[points_list[0].id, points_list[1].id], [points_list[1].id, points_list[2].id], [points_list[2].id, points_list[0].id]])

    lines = np.array([
        [0, 1],
        [1, 2],
        [2, 0],
        [1, 3],
        [3, 4],
        [4, 1],
        [3, 5],
        [5, 6],
        [6, 3],
        [5, 7],
        [7, 8],
        [8, 6],
        [7, 9],
        [9, 10],
        [10, 8],
        [9, 11],
        [11, 12],
        [12, 10],
        [2, 13],
        [13, 0],
        [4, 14],
        [14, 2],
        [6, 15],
        [15, 4],
        [8, 16],
        [16, 15],
    ])

    lines = lines+1

    lines2 = np.array([
        [11, 18],
        [18, 17],
        [13, 19],
        [19, 18],
        [15, 20],
        [20, 14],
        [16, 21],
        [21, 15],
        [17, 22],
        [22, 21],
        [18, 23],
        [23, 22],
        [19, 24],
        [24, 23],
        [20, 25],
        [25, 14],
        [21, 26],
        [26, 20],
        [22, 27], 
        [27, 26],
        [23, 28],
        [28, 27],
        [24, 29],
        [29, 28],
        [26, 30],
        [30, 25],
        [27, 31],
        [31, 30],
        [28, 32],
        [32, 31],
        [29, 33],
        [33, 32]
    ])

    lines = np.vstack([lines, lines2])
    
    curve_loop = np.array([
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9],
        [10, 11, 12, -8],
        [13, 14, 15, -11],
        [16, 17, 18, -14],
        [19, 20, -3],
        [21, 22, -2, -6],
        [23, 24, -5, -9],
        [25, 26, -23, -12],
        [27, 28, -25, -15],
        [29, 30, -27, -18],
        [31, 32, -19, -22],
        [33, 34, -21, -24],
        [35, 36, -33, -26],
        [37, 38, -35, -28],
        [39, 40, -37, -30],
        [41, 42, -32],
        [43, 44, -31, -34],
        [45, 46, -43, -36],
        [47, 48, -45, -38],
        [49, 50, -47, -40],
        [51, 52, -41, -44],
        [53, 54, -51, -46],
        [55, 56, -53, -48],
        [57, 58, -55, -50]
    ], dtype='O')

    elements = np.array([
        [1, 2, 3],
        [2, 4, 5],
        [4, 6, 7],
        [6, 8, 9, 7],
        [8, 10, 11, 9],
        [10, 12, 13, 11],
        [1, 3, 14],
        [2, 5, 15, 3],
        [5, 4, 7, 16],
        [7, 9, 17, 16],
        [9, 11, 18, 17],
        [11, 13, 19, 18],
        [14, 3, 15, 20],
        [15, 5, 16, 21],
        [21, 16, 17, 22],
        [22, 17, 18, 23],
        [23, 18, 19, 24],
        [25, 14, 20],
        [20, 15, 21, 26],
        [21, 22, 27, 26],
        [22, 23, 28, 27],
        [28, 23, 24, 29],
        [25, 20, 26, 30],
        [26, 27, 31, 30]
    ], dtype='O')

    text_points = write_points(points_list, lengths, points_ids_of_lenghts)
    text_lines = write_lines(lines)
    text_curves = write_curve_loop(curve_loop)

    write_file_geo(filename, text_points, text_lines, text_curves)

