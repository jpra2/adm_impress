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

    newpoints3 = newpoints1[0:4].copy()
    newpoints3[:, 1] = newpoints3[:, 1] + Ly/6
    points = np.vstack([points, newpoints3])

    newpoints4 = newpoints3.copy()
    newpoints4[:, 1] = newpoints4[:, 1] + Ly/12
    newpoints4[:, 0] = newpoints4[:, 0] - Lx/12
    points = np.vstack([points, newpoints4])

    newpoints5 = newpoints3[0:3].copy()
    newpoints5[:, 1] = newpoints5[:, 1] + Lx/6
    points = np.vstack([points, newpoints5])

    newpoints6 = newpoints5.copy()
    newpoints6[:, 1] = newpoints6[:, 1] + Lx/12
    newpoints6[:, 0] = newpoints6[:, 0] - Lx/12
    points = np.vstack([points, newpoints6])

    newpoints7 = newpoints6[1:].copy()
    newpoints7[:, 1] = newpoints7[:, 1] + Ly/12
    newpoints7[:, 0] = newpoints7[:, 0] - Lx/12
    points = np.vstack([points, newpoints7])

    newpoints8 = newpoints7.copy()
    newpoints8[:, 1] = newpoints8[:, 1] + Ly/12
    newpoints8[:, 0] = newpoints8[:, 0] - Lx/12
    points = np.vstack([points, newpoints8])

    newpoints9 = newpoints8[1:].copy()
    newpoints9[:, 1] = newpoints9[:, 1] + Ly/12
    newpoints9[:, 0] = newpoints9[:, 0] - Lx/12
    points = np.vstack([points, newpoints9])

    newpoints10 = newpoints9.copy()
    newpoints10[:, 1] = newpoints10[:, 1] + Ly/12
    newpoints10[:, 0] = newpoints10[:, 0] - Lx/12
    points = np.vstack([points, newpoints10])

    newpoints11 = newpoints10.copy()
    newpoints11[:, 0] = newpoints11[:, 0] + Lx/6
    points = np.vstack([points, newpoints11])

    points = np.hstack([points, np.zeros((points.shape[0], 1))])

    p51 = np.array([points[47].copy()])
    p51[:, 0] = p51[:, 0] + Lx/6
    points = np.vstack([points, p51])

    newpoints12 = points[[44, 50]].copy()
    newpoints12[:, 0] = newpoints12[:, 0] + Lx/12
    newpoints12[:, 1] = newpoints12[:, 1] + Ly/12
    points = np.vstack([points, newpoints12])

    newpoints13 = points[[44, 50]].copy()
    newpoints13[:, 0] = newpoints13[:, 0] + Lx/6
    points = np.vstack([points, newpoints13])

    newpoints14 = points[[42, 51, 52]].copy()
    newpoints14[:, 0] = newpoints14[:, 0] + Lx/6
    points = np.vstack([points, newpoints14])

    newpoints15 = points[[39, 53, 54]].copy()
    newpoints15[:, 0] = newpoints15[:, 0] + Lx/6
    points = np.vstack([points, newpoints15])

    newpoints16 = points[[36, 55, 56]].copy()
    newpoints16[:, 0] = newpoints16[:, 0] + Lx/6
    points = np.vstack([points, newpoints16])

    newpoints17 = points[[32, 58, 59]].copy()
    newpoints17[:, 0] = newpoints17[:, 0] + Lx/6
    points = np.vstack([points, newpoints17])

    newpoints18 = points[[28, 61, 62]].copy()
    newpoints18[:, 0] = newpoints18[:, 0] + Lx/6
    points = np.vstack([points, newpoints18])

    newpoints19 = points[[23, 64, 65]].copy()
    newpoints19[:, 0] = newpoints19[:, 0] + Lx/6
    points = np.vstack([points, newpoints19])

    newpoints20 = points[[18, 67, 68]].copy()
    newpoints20[:, 0] = newpoints20[:, 0] + Lx/6
    points = np.vstack([points, newpoints20])

    newpoints21 = points[[57]].copy()
    newpoints21[:, 0] = newpoints21[:, 0] + Lx/6
    points = np.vstack([points, newpoints21])

    newpoints22 = points[[60]].copy()
    newpoints22[:, 0] = newpoints22[:, 0] + Lx/6
    points = np.vstack([points, newpoints22])

    newpoints22 = points[[63]].copy()
    newpoints22[:, 0] = newpoints22[:, 0] + Lx/6
    points = np.vstack([points, newpoints22])

    newpoints23 = points[[66]].copy()
    newpoints23[:, 0] = newpoints23[:, 0] + Lx/6
    points = np.vstack([points, newpoints23])

    newpoints24 = points[[69]].copy()
    newpoints24[:, 0] = newpoints24[:, 0] + Lx/6
    points = np.vstack([points, newpoints24])

    newpoints25 = points[[76, 77, 78]].copy()
    newpoints25[:, 0] = newpoints25[:, 0] + Lx/6
    points = np.vstack([points, newpoints25])

    newpoints26 = points[[81]].copy()
    newpoints26[:, 0] = newpoints26[:, 0] + Lx/6
    points = np.vstack([points, newpoints26])

    lengths = [0.02, 0.01]

    points_list = np.array([Point(coord[0], coord[1], coord[2]) for coord in points])
    all_points_ids = np.array([point.id for point in points_list])
    points_ids_lenght2 = np.array([0, 2, 14, 20, 21, 22, 23, 64, 65, 66, 78, 82, 84,
                                   13, 19, 25, 26, 27, 28, 61, 62, 63, 77, 81,
                                   1, 4, 15, 16, 17, 18, 67, 68, 69, 79, 83,
                                   3, 24, 33, 40, 45, 48, 49, 52, 57, 80])
    points_ids_lenght1 = np.setdiff1d(all_points_ids, points_ids_lenght2)

    points_ids_of_lenghts = np.array([
        points_ids_lenght1,
        points_ids_lenght2
    ], dtype='O')



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
        [33, 32],
        [30, 34],
        [34, 25],
        [31, 35],
        [35, 34],
        [32, 36],
        [36, 35],
        [33, 37],
        [37, 36],
        [36, 39],
        [39, 38],
        [38, 35],
        [37, 40],
        [40, 39],
        [38, 34],
        [38, 41],
        [41, 34],
        [39, 42],
        [42, 41],
        [40, 43],
        [43, 42],
        [43, 45],
        [45, 44],
        [44, 42],
        [44, 41],
        [45, 47],
        [47, 46],
        [46, 44],
        [46, 41],
        [47, 48],
        [48, 46],
        [48, 49],
        [49, 46],
        [47, 50],
        [50, 48],
        [50, 49],
        [45, 51],
        [51, 50],
        [43, 52],
        [52, 51],
        [52, 53],
        [53, 51],
        [53, 50],
        [40, 54],
        [54, 52],
        [54, 55],
        [55, 53],
        [37, 56],
        [56, 54],
        [56, 57],
        [57, 55],
        [57, 58],
        [58, 55],
        [58, 53],
        [33, 59],
        [59, 56],
        [59, 60],
        [60, 57],
        [60, 61],
        [61, 58],
        [29, 62],
        [62, 59],
        [62, 63],
        [63, 60],
        [63, 64],
        [64, 61],
        [24, 65],
        [65, 62],
        [65, 66],
        [66, 63],
        [66, 67],
        [67, 64],
        [19, 68],
        [68, 65],
        [68, 69],
        [69, 66],
        [69, 70],
        [70, 67],
        [13, 71],
        [71, 68],
        [71, 72],
        [72, 69],
        [72, 73],
        [73, 70],
        [12, 74],
        [74, 71],
        [74, 75],
        [75, 72],
        [75, 76],
        [76, 73],
        [61, 77],
        [77, 58],
        [64, 78],
        [78, 77],
        [67, 79],
        [79, 78],
        [70, 80],
        [80, 79],
        [73, 81],
        [81, 80],
        [76, 81],
        [78, 82],
        [82, 77],
        [79, 83],
        [83, 82],
        [80, 84],
        [84, 83],
        [81, 84],
        [83, 85],
        [85, 82],
        [84, 85]
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
        [57, 58, -55, -50],
        [59, 60, -52],
        [61, 62, -59, -54],
        [63, 64, -61, -56],
        [65, 66, -63, -58],
        [67, 68, 69, -64],
        [70, 71, -67, -66],
        [-69, 72, -62],
        [73, 74, -72],
        [75, 76, -73, -68],
        [77, 78, -75, -71],
        [79, 80, 81, -78],
        [-81, 82, -76],
        [83, 84, 85, -80],
        [-85, 86, -82],
        [87, 88, -84],
        [89, 90, -88],
        [91, 92, -87],
        [93, -89, -92],
        [94, 95, -91, -83],
        [96, 97, -94, -79],
        [98, 99, -97],
        [100, -95, -99],
        [101, 102, -96, -77],
        [103, 104, -98, -102],
        [105, 106, -101, -70],
        [107, 108, -103, -106],
        [109, 110, -108],
        [111, -104, -110],
        [112, 113, -105, -65],
        [114, 115, -107, -113],
        [116, 117, -109, -115],
        [118, 119, -112, -57],
        [120, 121, -114, -119],
        [122, 123, -116, -121],
        [124, 125, -118, -49],
        [126, 127, -120, -125],
        [128, 129, -122, -127],
        [130, 131, -124, -39],
        [132, 133, -126, -131],
        [134, 135, -128, -133],
        [136, 137, -130, -29],
        [138, 139, -132, -137],
        [140, 141, -134, -139],
        [142, 143, -136, -17],
        [144, 145, -138, -143],
        [146, 147, -140, -145],
        [148, 149, -117],
        [150, 151, -148, -123],
        [152, 153, -150, -129],
        [154, 155, -152, -135],
        [156, 157, -154, -141],
        [158, -156, -147],
        [159, 160, -151],
        [161, 162, -159, -153],
        [163, 164, -161, -155],
        [165, -163, -157],
        [166, 167, -162],
        [168, -166, -164]
        
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
        [26, 27, 31, 30],
        [34, 25, 30],
        [30, 31, 35, 34],
        [31, 32, 36, 35],
        [32, 33, 37, 36],
        [35, 36, 39, 38],
        [36, 37, 40, 39],
        [34, 35, 38],
        [34, 38, 41],
        [38, 39, 42, 41],
        [39, 40, 43, 42],
        [42, 43, 45, 44],
        [41, 42, 44],
        [44, 45, 47, 46],
        [41, 44, 46],
        [46, 47, 48],
        [46, 48, 49],
        [48, 47, 50],
        [50, 49, 48],
        [47, 45, 51, 50],
        [45, 43, 52, 51],
        [51, 52, 53],
        [50, 51, 53],
        [52, 43, 40, 54],
        [53, 52, 54, 55],
        [54, 40, 37, 56],
        [55, 54, 56, 57],
        [55, 57, 58],
        [53, 55, 58],
        [56, 37, 33, 59],
        [57, 56, 59, 60],
        [58, 57, 60, 61],
        [59, 33, 29, 62],
        [60, 59, 62, 63],
        [61, 60, 63, 64],
        [65, 62, 29, 24],
        [64, 63, 66, 67],
        [65, 24, 19, 68],
        [66, 65, 68, 69],
        [67, 66, 69, 70],
        [68, 19, 13, 71],
        [69, 68, 71, 72],
        [70, 69, 72, 73],
        [71, 13, 12, 74],
        [72, 71, 74, 75],
        [73, 72, 75, 76],
        [58, 61, 77],
        [77, 61, 64, 78],
        [78, 64, 67, 79],
        [79, 67, 70, 80],
        [80, 70, 73, 81],
        [81, 73, 76],
        [77, 78, 82],
        [82, 78, 79, 83],
        [83, 79, 80, 84],
        [84, 80, 81],
        [85, 82, 83],
        [85, 83, 84]
    ], dtype='O')

    text_points = write_points(points_list, lengths, points_ids_of_lenghts)
    text_lines = write_lines(lines)
    text_curves = write_curve_loop(curve_loop)

    write_file_geo(filename, text_points, text_lines, text_curves)

