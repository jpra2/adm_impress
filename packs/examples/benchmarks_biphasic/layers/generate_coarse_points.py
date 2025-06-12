import numpy as np
from pymoab import core, types, rng, topo_util
from typing import Sequence
from packs.utils.calculate_face_properties import sort_radial_sweep, sort_vertices_by_zdirection_xy_plane
import os
from packs.mpfa_methods.mesh_preprocess import preprocess_mesh
import copy

def get_points_of_line(x1, y1, x2, y2, x) -> float:
    """return the y point of a line"""

    a = (y2 - y1)/(x2 - x1)
    b = y1 - a*x1

    y = a*x + b
    return y

class Line:
    
    def __init__(self, x1, y1, x2, y2):
        self.a = (y2 - y1)/(x2 - x1)
        self.b = y1 - self.a*x1
    
    def get_y(self, x):
        y = self.a*x + self.b
        return y
    
    def get_x(self, y):
        x = (y - self.b)/self.a
        return x

class Point:
    sequence = []
    def __init__(self, i: float, j: float, k: float):
        self.coordinate = np.array([i, j, k], dtype=np.float64)
        self.idp = len(Point.sequence)
        Point.sequence.append(self)
    
    def getx(self):
        return self.coord.copy()[0]
    
    def gety(self):
        return self.coord.copy()[1]
    
    def getz(self):
        return self.coord.copy()[2]
    
    @property
    def id(self):
        return copy.deepcopy(self.idp)

    @property
    def coord(self):
        return self.coordinate.copy()

    def show_name(self):
        return 'Point_' + str(self.id)

    def __str__(self):
        return self.show_name()


class Polygon:
    sequence = []
    def __init__(self, points: Sequence[Point]):
        self.points = points
        if len(self.points) == 3:
            self.type = types.MBTRI
        elif len(self.points) == 4:
            self.type = types.MBQUAD
        else:
            raise TypeError
        
        self.idp = len(Polygon.sequence)
        Polygon.sequence.append(self)

    @property
    def id(self):
        return copy.deepcopy(self.idp)

    @property
    def points_coords(self):
        p_coords = np.array(([point.coord for point in self.points]))
        return p_coords
    
    @property
    def points_ids(self):
        return np.array([point.id for point in self.points])
    
    @property
    def npoints(self):
        return len(self.points)

    
    def sort_by_xy_plane(self):
        p_coords = self.points_coords
        indexes = sort_radial_sweep(p_coords, np.arange(len(p_coords)))
        p_coords = p_coords[indexes]
        self.points = self.points[indexes]
        indexes2 = sort_vertices_by_zdirection_xy_plane(p_coords)
        self.points = self.points[indexes2]
    
    def point_intersection(self, polygon):
        polygon: Polygon
        if isinstance(polygon, Polygon):
            my_ids = self.points_ids
            another_ids = polygon.points_ids
            intersect = np.intersect1d(my_ids, another_ids)
            return intersect
        else:
            raise TypeError
    
    @property
    def centroid(self):
        p_coords = self.points_coords
        centroid = np.mean(p_coords, axis=0)
        return centroid








class MoabMesh:
    def __init__(self):
        self.mb = core.Core()
        self.root_set = self.mb.get_root_set()
        self.mtu = topo_util.MeshTopoUtil(self.mb)
        self.tri_type = types.MBTRI
        self.quad_type = types.MBQUAD
        self.poly_type = types.MBPOLYGON
    
    def create_vertexes(self, coords: np.ndarray):
        self.verts = self.mb.create_vertices(coords.flatten())
    
    def create_quads_elements(self, quads):
        self.mb.create_elements(self.quad_type, quads)
    
    def create_quad_element(self, quad):
        self.mb.create_element(self.quad_type, quad)
    
    def create_tri_elements(self, tris):
        self.mb.create_elements(self.tri_type, tris)
    
    def create_tri_element(self, tri):
        self.mb.create_element(self.tri_type, tri)
    
    def create_polygon(self, polygon):
        self.mb.create_element(self.poly_type, polygon)
    
    def export_mesh(self, mesh_name):
        mesh_name_export = os.path.join('mesh', mesh_name)
        self.mb.write_file(mesh_name_export)
    



def run():
    ## reta 1
    r1 = Line(0, 0.5, 0.5, 0)
    ## reta 2
    r2 = Line(0, 1, 1, 0)
    ## reta 3
    r3 = Line(0.5, 1, 1, 0.5)

    ## reta 4
    r4 = Line(0, 0.75, 0.75, 0)

    ## reta 5
    r5 = Line(0.25, 1, 1, 0.25)


    x0 = 0
    y0 = 0
    p0 = [x0,y0]

    x1 = 0.25
    y1 = r1.get_y(0.25)
    p1 = [x1, y1]

    x2 = (0.25 + 0.5)/2
    y2 = r1.get_y(x2)
    p2 = [x2, y2]

    x3 = 0.25
    y3 = 0.0
    p3 = [x3, y3]

    el1 = [p0, p1, p2, p3]

    x4 = 0.25/2
    y4 = r1.get_y(x4)
    p4 = [x4, y4]

    x5 = 0.0
    y5 = 0.25
    p5 = [x5, y5]

    x6 = 0.5
    y6 = 0.0
    p6 = [x6, y6]

    x7 = 0.0
    y7 = 0.5
    p7 = [x7, y7]

    y8 = 0.125
    x8 = r2.get_x(y8)
    p8 = [x8, y8]

    x9 = 1.0
    y9 = 0.0
    p9 = [x9, y9]

    x10 = 0.625
    y10 = 0.125
    p10 = [x10, y10]

    x11 = 0.75
    y11 = 0.0
    p11 = [x11, y11]

    y12 = 0.25
    x12 = r2.get_x(y12)
    p12 = [x12, y12]

    y13 = 0.25
    x13 = 0.5
    p13 = [x13, y13]

    y14 = 0.375
    x14 = r2.get_x(y14)
    p14 = [x14, y14]

    x15 = (x4 + x14)/2
    y15 = y4
    p15 = [x15, y15]

    x16 = x4
    y16 = y2
    p16 = [x16, y16]

    p17 = [x13, y7]

    p18 = [x1, y7]

    p19 = [r2.get_x(0.625), 0.625]
    p20 = [p4[0], p19[1]]
    p21 = [r2.get_x(0.75), 0.75]

    p22 = [0, 0.75]
    p23 = [r2.get_x(0.875), 0.875]
    p24 = [0, 1]
    p25 = [0.25, 1]
    p26 = [0.375, 0.875]
    p27 = [0.5, 1.0]
    p28 = [0.5, 0.75]
    p29 = [0.625, 0.625]
    p30 = [0.625, 0.875]
    p31 = [0.75, 0.5]
    p32 = [0.75, 0.75]
    p33 = [0.875, 0.375]
    p34 = [0.875, 0.625]
    p35 = [1, 0.25]
    p36 = [1, 0.5]
    p37 = [1, 0.75]
    p38 = [0.875, 0.875]
    p39 = [0.75, 1]
    p40 = [1, 1]

    points = [p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16,
              p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29,
              p30, p31, p32, p33, p34, p35, p36, p37, p38, p39, p40]
    
    points = np.array(points)
    points = np.hstack([points, np.zeros((points.shape[0], 1))])

    elements = [
        [0, 3, 16],
        [0, 16, 5],
        [3, 6, 2],
        [16, 3, 2, 1],
        [5, 16, 1, 4],
        [5, 4, 7],
        [2, 6, 11, 10],
        [10, 8, 9, 11],
        [2, 10, 13, 1],
        [10, 8, 12, 13],
        [1, 13, 15, 4],
        [13, 12, 14, 15],
        [4, 15, 18, 7],
        [15, 14, 17, 18],
        [7, 18, 20],
        [18, 17, 19, 20],
        [7, 20, 22],
        [20, 19, 21, 22],
        [22, 21, 23],
        [22, 23, 24],
        [23, 25, 24],
        [23, 21, 25],
        [21, 19, 26, 25],
        [26, 27, 25],
        [19, 17, 28, 26],
        [26, 28, 27],
        [17, 14, 29, 28],
        [28, 29, 30, 27],
        [14, 12, 31, 29],
        [29, 31, 32, 30],
        [12, 8, 33, 31],
        [31, 33, 34, 32],
        [8, 9, 35, 33],
        [33, 35, 36, 34],
        [34, 36, 37],
        [32, 34, 37, 38],
        [30, 32, 38, 39],
        [38, 37, 40],
        [39, 38, 40],
        [27, 30, 39]
    ]

    points_list = np.array([Point(coord[0], coord[1], coord[2]) for coord in points])
    polygons = np.array([Polygon(points_list[indexes]) for indexes in elements])
    for i, polygon in enumerate(polygons):
        polygon: Polygon
        polygon.sort_by_xy_plane()
    
    polygons_n_points = np.array([poly.npoints for poly in polygons])

    polygons3 = polygons[polygons_n_points==3]
    polygons4 = polygons[polygons_n_points==4]
    
    moab = MoabMesh()
    moab.create_vertexes(points)

    quads = [moab.verts[polygon.points_ids] for polygon in polygons4]
    tris = [moab.verts[polygon.points_ids] for polygon in polygons3]

    moab.create_tri_elements(tris)
    # moab.create_quads_elements(quads)
    print(len(tris))
    import pdb; pdb.set_trace()

    list_quads = [2, 3, 5, 6, 7, 8, 9, 10, 11, 17, 19, 20, 21]

    for i, quad in enumerate(quads):
        moab.create_quad_element(quad)
        moab.export_mesh('coarse_test.vtk')
        moab.export_mesh('coarse_test.msh')
        print(i)
        import pdb; pdb.set_trace()


    mesh = preprocess_mesh('coarse_test.msh', 'coarse_test')
    import pdb; pdb.set_trace()

    

    # for i, quad in enumerate(quads):
    #     moab.create_quad_element(quad)
    #     moab.export_mesh('coarse_test.vtk')
        

    import pdb; pdb.set_trace()

        

    # moab.create_quads_elements(quads)
    

    # import pdb; pdb.set_trace()

def run2():
    import meshio

    # two triangles and one quad
    points = [
        [0.0, 0.0],
        [1.0, 0.0],
        [0.0, 1.0],
        [1.0, 1.0],
        [2.0, 0.0],
        [2.0, 1.0],
    ]
    cells = [
        ("triangle", [[0, 1, 2], [1, 3, 2]]),
        ("quad", [[1, 4, 5, 3]]),
    ]

    mesh = meshio.Mesh(
        points,
        cells,
        # Optionally provide extra data on points, cells, etc.
        point_data={"T": [0.3, -1.2, 0.5, 0.7, 0.0, -3.0]},
        # Each item in cell data must match the cells array
        cell_data={"a": [[0.1, 0.2], [0.4]]},
    )
    # mesh.write(
    #     "foo.vtk",  # str, os.PathLike, or buffer/open file
    #     # file_format="vtk",  # optional if first argument is a path; inferred from extension
    # )

    meshio.write_points_cells("foo.vtk", points, cells)

    ##
    """
    $MeshFormat
    2.2 0 8
    $EndMeshFormat
    $Nodes
    Numero total de nos
    nome do no, coordenadas x, y, z 
    $EndNodes
    $Elements
    Numero total de elementos
    """



