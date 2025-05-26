import numpy as np
import os
from packs import defpaths

def pstr(npoint, lelement, x, y):

    val = 'Point(' + str(npoint) + ') = {' + str(x) + ',' + str(y) + ',0,' + str(lelement) + '};'
    return val

def bspline_str(points_ids, ncurve):

    points_str = [str(i) for i in points_ids]
    points_str = ','.join(points_str)
    # points_str = points_str + ',' + str(points_ids[0])
    points_str = '{' + points_str + '};'
    bspline = 'BSpline(' + str(ncurve) + ') = ' + points_str
    return bspline

def point_in(point_id, plane_surface):
    # Point{21} In Surface {1};
    v1 = 'Point{' + str(point_id) + '} In Surface {' + str(plane_surface) + '};'
    return v1

def chueh_sin_geo_coarse():

    cl1 = 0.25
    cl2 = 0.3

    x = np.linspace(0, 1, 51)
    nx = x.shape[0]
    y = np.round(0.5 + 0.1*np.sin(10*x), 6)
    y2 = y + 0.12
    y3 = y - 0.12

    text = []    

    for i in range(nx):
        v1 = pstr(i+1, cl1, x[i], y2[i])
        v2 = pstr(i+nx+1, cl1, x[i], y3[i])
        text.append(v1)
        text.append(v2)
    
    bs1 = bspline_str(np.arange(1, nx+1), 1)
    bs2 = bspline_str(np.arange(nx+1, 2*nx+1), 2)
    text.append(bs1)
    text.append(bs2)

    total_points = 2*nx
    tp = total_points

    p0 = pstr(total_points+1, cl2, 0, 0)
    p1 = pstr(total_points+2, cl2, 1, 0)
    p2 = pstr(total_points+3, cl2, 1, 1)
    p3 = pstr(total_points+4, cl2, 0, 1)
    lopints2 = [p0, p1, p2, p3]
    for i in lopints2:
        text.append(i)
    
    nlines = 2

    l1 = 'Line(' + str(nlines+1) + ')={' + str(tp+1) + ',' + str(tp+2) + '};'
    l2 = 'Line(' + str(nlines+2) + ')={' + str(tp+2) + ',' + str(2*nx) + '};'
    l3 = 'Line(' + str(nlines+3) + ')={' + str(2*nx) + ',' + str(nx) + '};'
    l4 = 'Line(' + str(nlines+4) + ')={' + str(nx) + ',' + str(tp+3) + '};'
    l5 = 'Line(' + str(nlines+5) + ')={' + str(tp+3) + ',' + str(tp+4) + '};'
    l6 = 'Line(' + str(nlines+6) + ')={' + str(tp+4) + ',' + str(1) + '};'
    l7 = 'Line(' + str(nlines+7) + ')={' + str(1) + ',' + str(nx+1) + '};'
    l8 = 'Line(' + str(nlines+8) + ')={' + str(nx+1) + ',' + str(tp+1) + '};'
    l9 = '\nTransfinite Line(' + str(1) + ')= ' + str(nx-25) + ' Using Progression 1;'
    l10 = 'Transfinite Line(' + str(2) + ')= ' + str(nx-25) + ' Using Progression 1;\n'

    nlineloop = 1
    
    l11 = 'Line Loop(1) = {1, 6, 7, 8};' 
    l12 = 'Line Loop(2) = {9, 2, 5, -1};' 
    l13 = 'Line Loop(3) = {10, 3, 4, -2};'
    l14 = 'Plane Surface(1)={1};' 
    l15 = 'Plane Surface(2)={2};' 
    l16 = 'Plane Surface(3)={3};' 

    # l11 = 'Line Loop(1) = {3, 4, 5, 6, 7, 8, 9, 10};'
    # l12 = 'Line Loop(2) = {9, 2, 5, -1};'
    # l14 = 'Plane Surface(1) = {1, -2};'
    # l15 = 'Plane Surface(2) = {2};'


    l17 = 'Mesh.Algorithm =2;'
    l18 = 'Physical Line(201) = {8, 9, 10};'
    l19 = 'Physical Line(202) = {4, 5, 6};'
    l20 = 'Physical Line(203) = {3, 7};'
    l21 = 'Recombine Surface{1, 2, 3};'


    # llines = [l1, l2, l3, l4, l5, l6, l7, l8, l11, l12, l13, l14, l15, l16, l17,
    #           l18, l19, l20]
    
    llines = [l1, l2, l3, l4, l5, l6, l7, l8, l11, l12, l13, l14, l15, l16, l17,
              l18, l19, l20]

    llines = [l1, l2, l3, l4, l5, l6, l7, l8, l11, l12, l13, l14, l15, l16, l17,
              l18, l19, l20, l21]
    
    # llines = [l1, l2, l3, l4, l5, l6, l7, l8, l11, l12, l14, l15, l17,
    #           l18, l19]



    for i in llines:
        text.append(i)

    # for pid in np.arange(1, 2*nx+1):
    #     v1 = point_in(pid, 2)
    #     text.append(v1)
    
    # for pid in np.arange(1, nx+1):
    #     v1 = point_in(pid, 1)
    #     text.append(v1)
    
    # for pid in np.arange(nx+1, 2*nx+1):
    #     v1 = point_in(pid, 3)
    #     text.append(v1)
    
    file_path = os.path.join(defpaths.mesh, defpaths.biphasic_folder, 'sin_coarse.geo')

    with open(file_path, 'w') as f:
        for line in text:
            lline = line + '\n'
            f.write(lline)
    


    
    
    import pdb; pdb.set_trace()

def chueh_sin_geo_coarse2():

    cl1 = 0.15
    cl2 = 0.2

    x = np.linspace(0, 1, 51)
    nx = x.shape[0]
    y = np.round(0.5 + 0.1*np.sin(10*x), 6)
    y2 = y + 0.12
    y3 = y - 0.12

    text = []    

    for i in range(nx):
        v1 = pstr(i+1, cl1, x[i], y2[i])
        v2 = pstr(i+nx+1, cl1, x[i], y3[i])
        text.append(v1)
        text.append(v2)
    
    bs1 = bspline_str(np.arange(1, nx+1), 1)
    bs2 = bspline_str(np.arange(nx+1, 2*nx+1), 2)
    text.append(bs1)
    text.append(bs2)

    total_points = 2*nx
    tp = total_points

    p0 = pstr(total_points+1, cl2, 0, 0)
    p1 = pstr(total_points+2, cl2, 1, 0)
    p2 = pstr(total_points+3, cl2, 1, 1)
    p3 = pstr(total_points+4, cl2, 0, 1)
    lopints2 = [p0, p1, p2, p3]
    for i in lopints2:
        text.append(i)
    
    nlines = 2

    l1 = 'Line(' + str(nlines+1) + ')={' + str(tp+1) + ',' + str(tp+2) + '};'
    l2 = 'Line(' + str(nlines+2) + ')={' + str(tp+2) + ',' + str(2*nx) + '};'
    l3 = 'Line(' + str(nlines+3) + ')={' + str(2*nx) + ',' + str(nx) + '};'
    l4 = 'Line(' + str(nlines+4) + ')={' + str(nx) + ',' + str(tp+3) + '};'
    l5 = 'Line(' + str(nlines+5) + ')={' + str(tp+3) + ',' + str(tp+4) + '};'
    l6 = 'Line(' + str(nlines+6) + ')={' + str(tp+4) + ',' + str(1) + '};'
    l7 = 'Line(' + str(nlines+7) + ')={' + str(1) + ',' + str(nx+1) + '};'
    l8 = 'Line(' + str(nlines+8) + ')={' + str(nx+1) + ',' + str(tp+1) + '};'
    l9 = '\nTransfinite Line(' + str(1) + ')= ' + str(nx-25) + ' Using Progression 1;'
    l10 = 'Transfinite Line(' + str(2) + ')= ' + str(nx-25) + ' Using Progression 1;\n'

    nlineloop = 1
    
    l11 = 'Line Loop(1) = {1, 6, 7, 8};' 
    l12 = 'Line Loop(2) = {9, 2, 5, -1};' 
    l13 = 'Line Loop(3) = {10, 3, 4, -2};'
    l14 = 'Plane Surface(1)={1};' 
    l15 = 'Plane Surface(2)={2};' 
    l16 = 'Plane Surface(3)={3};' 

    # l11 = 'Line Loop(1) = {3, 4, 5, 6, 7, 8, 9, 10};'
    # l12 = 'Line Loop(2) = {9, 2, 5, -1};'
    # l14 = 'Plane Surface(1) = {1, -2};'
    # l15 = 'Plane Surface(2) = {2};'


    l17 = 'Mesh.Algorithm =2;'
    l18 = 'Physical Line(201) = {8, 9, 10};'
    l19 = 'Physical Line(202) = {4, 5, 6};'
    l20 = 'Physical Line(203) = {3, 7};'
    l21 = 'Recombine Surface{1, 2, 3};'


    # llines = [l1, l2, l3, l4, l5, l6, l7, l8, l11, l12, l13, l14, l15, l16, l17,
    #           l18, l19, l20]
    
    llines = [l1, l2, l3, l4, l5, l6, l7, l8, l11, l12, l13, l14, l15, l16, l17,
              l18, l19, l20]

    llines = [l1, l2, l3, l4, l5, l6, l7, l8, l11, l12, l13, l14, l15, l16, l17,
              l18, l19, l20, l21]
    
    # llines = [l1, l2, l3, l4, l5, l6, l7, l8, l11, l12, l14, l15, l17,
    #           l18, l19]



    for i in llines:
        text.append(i)

    # for pid in np.arange(1, 2*nx+1):
    #     v1 = point_in(pid, 2)
    #     text.append(v1)
    
    # for pid in np.arange(1, nx+1):
    #     v1 = point_in(pid, 1)
    #     text.append(v1)
    
    # for pid in np.arange(nx+1, 2*nx+1):
    #     v1 = point_in(pid, 3)
    #     text.append(v1)
    
    file_path = os.path.join(defpaths.mesh, defpaths.biphasic_folder, 'sin_coarse2.geo')

    with open(file_path, 'w') as f:
        for line in text:
            lline = line + '\n'
            f.write(lline)
    


    
    
    import pdb; pdb.set_trace()

def chueh_sin_geo_coarse3():

    cl1 = 0.15
    cl2 = 0.15

    x = np.linspace(0, 1, 51)
    nx = x.shape[0]
    y = np.round(0.5 + 0.1*np.sin(10*x), 6)
    y2 = y + 0.21
    y3 = y - 0.21

    text = []    

    for i in range(nx):
        v1 = pstr(i+1, cl1, x[i], y2[i])
        v2 = pstr(i+nx+1, cl1, x[i], y3[i])
        text.append(v1)
        text.append(v2)
    
    bs1 = bspline_str(np.arange(1, nx+1), 1)
    bs2 = bspline_str(np.arange(nx+1, 2*nx+1), 2)
    text.append(bs1)
    text.append(bs2)

    total_points = 2*nx
    tp = total_points

    p0 = pstr(total_points+1, cl2, 0, 0)
    p1 = pstr(total_points+2, cl2, 1, 0)
    p2 = pstr(total_points+3, cl2, 1, 1)
    p3 = pstr(total_points+4, cl2, 0, 1)
    lopints2 = [p0, p1, p2, p3]
    for i in lopints2:
        text.append(i)
    
    nlines = 2

    l1 = 'Line(' + str(nlines+1) + ')={' + str(tp+1) + ',' + str(tp+2) + '};'
    l2 = 'Line(' + str(nlines+2) + ')={' + str(tp+2) + ',' + str(2*nx) + '};'
    l3 = 'Line(' + str(nlines+3) + ')={' + str(2*nx) + ',' + str(nx) + '};'
    l4 = 'Line(' + str(nlines+4) + ')={' + str(nx) + ',' + str(tp+3) + '};'
    l5 = 'Line(' + str(nlines+5) + ')={' + str(tp+3) + ',' + str(tp+4) + '};'
    l6 = 'Line(' + str(nlines+6) + ')={' + str(tp+4) + ',' + str(1) + '};'
    l7 = 'Line(' + str(nlines+7) + ')={' + str(1) + ',' + str(nx+1) + '};'
    l8 = 'Line(' + str(nlines+8) + ')={' + str(nx+1) + ',' + str(tp+1) + '};'
    l9 = '\nTransfinite Line(' + str(1) + ')= ' + str(nx-25) + ' Using Progression 1;'
    l10 = 'Transfinite Line(' + str(2) + ')= ' + str(nx-25) + ' Using Progression 1;\n'

    nlineloop = 1
    
    l11 = 'Line Loop(1) = {1, 6, 7, 8};' 
    l12 = 'Line Loop(2) = {9, 2, 5, -1};' 
    l13 = 'Line Loop(3) = {10, 3, 4, -2};'
    l14 = 'Plane Surface(1)={1};' 
    l15 = 'Plane Surface(2)={2};' 
    l16 = 'Plane Surface(3)={3};' 

    # l11 = 'Line Loop(1) = {3, 4, 5, 6, 7, 8, 9, 10};'
    # l12 = 'Line Loop(2) = {9, 2, 5, -1};'
    # l14 = 'Plane Surface(1) = {1, -2};'
    # l15 = 'Plane Surface(2) = {2};'


    l17 = 'Mesh.Algorithm =2;'
    l18 = 'Physical Line(201) = {8, 9, 10};'
    l19 = 'Physical Line(202) = {4, 5, 6};'
    l20 = 'Physical Line(203) = {3, 7};'
    l21 = 'Recombine Surface{1, 2, 3};'


    # llines = [l1, l2, l3, l4, l5, l6, l7, l8, l11, l12, l13, l14, l15, l16, l17,
    #           l18, l19, l20]
    
    llines = [l1, l2, l3, l4, l5, l6, l7, l8, l11, l12, l13, l14, l15, l16, l17,
              l18, l19, l20]

    llines = [l1, l2, l3, l4, l5, l6, l7, l8, l11, l12, l13, l14, l15, l16, l17,
              l18, l19, l20, l21]
    
    # llines = [l1, l2, l3, l4, l5, l6, l7, l8, l11, l12, l14, l15, l17,
    #           l18, l19]



    for i in llines:
        text.append(i)

    # for pid in np.arange(1, 2*nx+1):
    #     v1 = point_in(pid, 2)
    #     text.append(v1)
    
    # for pid in np.arange(1, nx+1):
    #     v1 = point_in(pid, 1)
    #     text.append(v1)
    
    # for pid in np.arange(nx+1, 2*nx+1):
    #     v1 = point_in(pid, 3)
    #     text.append(v1)
    
    file_path = os.path.join(defpaths.mesh, defpaths.biphasic_folder, 'sin_coarse3.geo')

    with open(file_path, 'w') as f:
        for line in text:
            lline = line + '\n'
            f.write(lline)
    


    
    
    import pdb; pdb.set_trace()


def chueh_sin_geo():

    cl1 = 0.015
    cl2 = 0.015

    x = np.linspace(0, 1, 51)
    nx = x.shape[0]
    y = np.round(0.5 + 0.1*np.sin(10*x), 6)
    y2 = y + 0.12
    y3 = y - 0.12

    text = []    

    # for i in range(nx):
    #     v1 = pstr(i+1, cl1, x[i], y2[i])
    #     v2 = pstr(i+nx+1, cl1, x[i], y3[i])
    #     text.append(v1)
    #     text.append(v2)
    
    # bs1 = bspline_str(np.arange(1, nx+1), 1)
    # bs2 = bspline_str(np.arange(nx+1, 2*nx+1), 2)
    # text.append(bs1)
    # text.append(bs2)

    total_points = 2*nx
    tp = total_points

    p0 = pstr(total_points+1, cl2, 0, 0)
    p1 = pstr(total_points+2, cl2, 1, 0)
    p2 = pstr(total_points+3, cl2, 1, 1)
    p3 = pstr(total_points+4, cl2, 0, 1)
    lopints2 = [p0, p1, p2, p3]
    for i in lopints2:
        text.append(i)
    
    nlines = 2

    l1 = 'Line(' + str(nlines+1) + ')={' + str(tp+1) + ',' + str(tp+2) + '};'
    l2 = 'Line(' + str(nlines+2) + ')={' + str(tp+2) + ',' + str(tp+3) + '};'
    l3 = 'Line(' + str(nlines+3) + ')={' + str(tp+3) + ',' + str(tp+4) + '};'
    l4 = 'Line(' + str(nlines+4) + ')={' + str(tp+4) + ',' + str(tp+1) + '};'
    # l5 = 'Line(' + str(nlines+5) + ')={' + str(tp+3) + ',' + str(tp+4) + '};'
    # l6 = 'Line(' + str(nlines+6) + ')={' + str(tp+4) + ',' + str(1) + '};'
    # l7 = 'Line(' + str(nlines+7) + ')={' + str(1) + ',' + str(nx+1) + '};'
    # l8 = 'Line(' + str(nlines+8) + ')={' + str(nx+1) + ',' + str(tp+1) + '};'
    l9 = '\nTransfinite Line(' + str(6) + ')= ' + str(91) + ' Using Progression 1;'
    l10 = 'Transfinite Line(' + str(4) + ')= ' + str(91) + ' Using Progression 1;'
    l12 = 'Transfinite Line(' + str(3) + ')= ' + str(61) + ' Using Progression 1;'
    l13 = 'Transfinite Line(' + str(5) + ')= ' + str(61) + ' Using Progression 1\n;'

    nlineloop = 1
    
    l11 = 'Line Loop(1) = {3, 4, 5, 6};' 
    # l12 = 'Line Loop(2) = {9, 2, 5, -1};' 
    # l13 = 'Line Loop(3) = {10, 3, 4, -2};'
    l14 = 'Plane Surface(1)={1};' 
    # l15 = 'Plane Surface(2)={2};' 
    # l16 = 'Plane Surface(3)={3};' 

    # # Point{21} In Surface {1};
    # l15 = 'Line{1} In Surface{1};'
    # l16 = 'Line{2} In Surface{1};'

    # l11 = 'Line Loop(1) = {3, 4, 5, 6, 7, 8, 9, 10};'
    # l12 = 'Line Loop(2) = {9, 2, 5, -1};'
    # l14 = 'Plane Surface(1) = {1, -2};'
    # l15 = 'Plane Surface(2) = {2};'


    l17 = 'Mesh.Algorithm =2;'
    l18 = 'Physical Line(201) = {6};'
    l19 = 'Physical Line(202) = {4};'
    l20 = 'Physical Line(203) = {3, 5};'


    # llines = [l1, l2, l3, l4, l5, l6, l7, l8, l11, l12, l13, l14, l15, l16, l17,
    #           l18, l19, l20]
    
    # llines = [l1, l2, l3, l4, l5, l6, l7, l8, l11, l12, l13, l14, l15, l16, l17,
    #           l18, l19, l20]
    
    # llines = [l1, l2, l3, l4, l5, l6, l7, l8, l11, l12, l14, l15, l17,
    #           l18, l19]

    llines = [l1, l2, l3, l4, l9, l10, l11, l12, l13, l14, l17, l18, l19, l20]



    for i in llines:
        text.append(i)

    for pid in np.arange(1, 2*nx+1):
        v1 = point_in(pid, 1)
        text.append(v1)
    
    # for pid in np.arange(1, nx+1):
    #     v1 = point_in(pid, 1)
    #     text.append(v1)
    
    # for pid in np.arange(nx+1, 2*nx+1):
    #     v1 = point_in(pid, 3)
    #     text.append(v1)
    
    file_path = os.path.join(defpaths.mesh, defpaths.biphasic_folder, 'sin.geo')

    with open(file_path, 'w') as f:
        for line in text:
            lline = line + '\n'
            f.write(lline)
    


    
    
    import pdb; pdb.set_trace()







