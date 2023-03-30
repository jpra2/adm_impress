import numpy as np
import math

def sort_radial_sweep_by_centroid(centroid, points, points_indices):
    cent = centroid
    vs = points
    indices = points_indices
    
    r0 = (vs[0] - cent)/np.linalg.norm(vs[0] - cent)
    nor = np.cross(vs[1] - cent, r0)
    nor = nor/np.linalg.norm(nor)
    
    # Pairs of (vertex index, angle to centroid)
    vpairs = []
    for vi, vpos in zip(indices, vs):
        # r1 = (vpos - cent).normalized()
        r1 = vpos - cent
        r1 = r1/np.linalg.norm(r1)
        dot = np.dot(r1,r0)
        angle = math.acos(max(min(dot, 1), -1))
        # angle *= 1 if nor.dot(r1.cross(r0)) >= 0 else -1
        dottest = np.dot(nor, np.cross(r1, r0))
        angle *= 1 if dottest >= 0 else -1        
        vpairs.append((vi, angle))
    
    # Sort by angle and return indices
    vpairs.sort(key=lambda v: v[1])
    # return np.array([vi for vi, angle in vpairs])
    return np.array([vi for vi, angle in vpairs])
    

def sort_radial_sweep(vs, indices):
    """
    Given a list of vertex positions (vs) and indices
    for verts making up a circular-ish planar polygon,
    returns the vertex indices in order around that poly.
    """
    # indices = np.arange(len(vs))
    assert len(vs) >= 3
    
    # Centroid of verts
    # cent = Vector()
    # for v in vs:
    #     cent += (1/len(vs)) * v
    
    cent = np.mean(vs, axis=0)

    
    # Normalized vector from centroid to first vertex
    # ASSUMES: vs[0] is not located at the centroid
    # r0 = (vs[0] - cent).normalized()
    r0 = (vs[0] - cent)/np.linalg.norm(vs[0] - cent)

    # Normal to plane of poly
    # ASSUMES: cent, vs[0], and vs[1] are not colinear
    # nor = (vs[1] - cent).cross(r0).normalized()
    nor = np.cross(vs[1] - cent, r0)
    nor = nor/np.linalg.norm(nor)

    # Pairs of (vertex index, angle to centroid)
    vpairs = []
    for vi, vpos in zip(indices, vs):
        # r1 = (vpos - cent).normalized()
        r1 = vpos - cent
        r1 = r1/np.linalg.norm(r1)
        dot = np.dot(r1,r0)
        angle = math.acos(max(min(dot, 1), -1))
        # angle *= 1 if nor.dot(r1.cross(r0)) >= 0 else -1
        dottest = np.dot(nor, np.cross(r1, r0))
        angle *= 1 if dottest >= 0 else -1        
        vpairs.append((vi, angle))
    
    # Sort by angle and return indices
    vpairs.sort(key=lambda v: v[1])
    # return np.array([vi for vi, angle in vpairs])
    return np.array([vi for vi, angle in vpairs])

def sort_vertices_by_zdirection_xy_plane(vs, zdirection=np.array([0, 0, 1])):
    cent = np.mean(vs, axis=0)
    index = np.arange(len(vs))
    v1 = vs[0] - cent
    v2 = vs[1] - cent
    cross_vector = np.cross(v1, v2)
    dot_prod = np.dot(cross_vector, zdirection)
    if dot_prod > 0:
        return index
    else:
        return np.flip(index)

def sort_vertices_by_zdirection_xy_plane_by_centroid(centroid, vs, zdirection=np.array([0, 0, 1])):
    cent = centroid
    index = np.arange(len(vs))
    v1 = vs[0] - cent
    v2 = vs[1] - cent
    cross_vector = np.cross(v1, v2)
    dot_prod = np.dot(cross_vector, zdirection)
    if dot_prod > 0:
        return index
    else:
        return np.flip(index)

def polygon_area(poly):
    """Calculate the area of 3d polygon from ordered vertices points (poly)

    Args:
        poly (_type_): oredered vertices of polygon. #shape (N, 3)

    Returns:
        area: polygon area
    """
    if isinstance(poly, list):
        poly = np.array(poly)
    #all edges
    edges = poly[1:] - poly[0:1]
    # row wise cross product
    cross_product = np.cross(edges[:-1],edges[1:], axis=1)
    #area of all triangles
    area = np.linalg.norm(cross_product, axis=1)/2
    return sum(area)

def get_unitary_normal_vector(vs):
    """return unitary normal vector

    Args:
        vs (_type_): points of polygon vertices
    """
    cent = np.mean(vs, axis=0)
    v1 = vs[0] - cent
    v2 = vs[1] - cent
    normal = np.cross(v1, v2)
    unitary = normal/np.linalg.norm(normal)
    return unitary
    
def sort_vertices_of_all_faces(faces, vertices_of_faces, vertices_centroids):
    """reordena os vertices em ordem circular,
    os vertices sao alterados no proprio array
    
    ###########
        obs
        esse codigo pode ser paralelizado
    ###########

    Args:
        faces (_type_): face ids
        vertices_of_faces (_type_): _description_
    """
    
    for face in faces:
        vf = vertices_of_faces[face]
        indices = np.arange(len(vf))
        cent_vertices = vertices_centroids[vf]
        new_indices = sort_radial_sweep(cent_vertices, indices)
        vertices_of_faces[face][:] = vf[new_indices]   
    
def define_normal_and_area(faces, vertices_of_faces, vertices_centroids):
    """return the area and unitary normal

    Args:
        faces (_type_): faces ids
        vertices_of_faces (_type_): vertices of faces
        volumes_adj_by_faces (_type_): volumes adjacencies
        volumes_centroids (_type_): volumes centroids
    """
    
    all_unitary_normals = np.zeros((len(faces), 3))
    all_areas = np.zeros(len(faces))
    
    sort_vertices_of_all_faces(faces, vertices_of_faces, vertices_centroids)
    
    for face in faces:
        vf = vertices_of_faces[face]
        cent_vertices = vertices_centroids[vf]
        area = polygon_area(cent_vertices)
        unitary_normal = get_unitary_normal_vector(cent_vertices)
        all_unitary_normals[face][:] = unitary_normal
        all_areas[face] = area
    
    return all_areas, all_unitary_normals

def correction_faces_vertices_order(unitary_normal_vector, nodes_of_faces, vector):
    test = unitary_normal_vector == vector
    test = np.all(test, axis=1)
    test = ~test
    nodes_of_faces[test] = np.flip(nodes_of_faces[test], axis=1)
    
def create_unitary_normal_edges_xy_plane(nodes_of_edges, centroids_of_nodes, faces_adj_by_edges, faces_centroids, bool_boundary_edges):
    
    vertices_of_edges = centroids_of_nodes[nodes_of_edges]
    edges_centroids = np.mean(vertices_of_edges, axis=1)
    vector_edges = vertices_of_edges[:, 1] - vertices_of_edges[:, 0]
    
    
    ## matriz rotação de pi/2 no plano xy
    if vector_edges.shape[1] == 3:
        R_matrix = np.array([[0, 1, 0],
                             [-1, 0, 0],
                             [0, 0, 0]])
    elif vector_edges.shape[1] == 2:
        R_matrix = np.array([[0, 1],
                             [-1, 0]])
    else:
        raise ValueError
    
    
    faces_direction_vector = faces_centroids[faces_adj_by_edges[:, 1]] - faces_centroids[faces_adj_by_edges[:, 0]]
    
    boundary_edges = bool_boundary_edges
    
    faces_direction_vector[boundary_edges] = edges_centroids[boundary_edges] - faces_centroids[faces_adj_by_edges[boundary_edges, 0]]
    
    normal_edges = np.matmul(vector_edges, R_matrix)
    norm = np.linalg.norm(normal_edges, axis=1)
    unitary_normal_edges = normal_edges/norm.reshape((norm.shape[0], 1))
    
    proj = np.diag(np.tensordot(unitary_normal_edges, faces_direction_vector, axes=((1), (1))))
    
    test = proj < 0
    unitary_normal_edges[test] = -1*unitary_normal_edges[test]
    
    return norm, unitary_normal_edges

def distance_from_point_to_line(point, line_point_1, line_point_2):
    
    norm = np.linalg.norm
    p1 = line_point_1
    p2 = line_point_2

    p3 = point
    distance = np.abs(norm(np.cross(p2-p1, p1-p3)))/norm(p2-p1)
    
    return distance

def create_face_to_edge_distances(faces_centroids, faces_adj_by_edges, nodes_of_edges, edges, nodes_centroids, bool_boundary_edges):
    
    face_to_edge_distance = np.zeros(faces_adj_by_edges.shape)
    
    bool_internal_edges = ~bool_boundary_edges
    internal_edges = edges[bool_internal_edges]
    
    for edge, faces in zip(internal_edges, faces_adj_by_edges[bool_internal_edges]):
        centroids_of_nodes_edge = nodes_centroids[nodes_of_edges[edge]]
        face_to_edge_distance[edge][:] = [
            distance_from_point_to_line(faces_centroids[faces[0]], centroids_of_nodes_edge[0], centroids_of_nodes_edge[1]),
            distance_from_point_to_line(faces_centroids[faces[1]], centroids_of_nodes_edge[0], centroids_of_nodes_edge[1])
        ]
    
    bedges = edges[bool_boundary_edges]
    
    for edge, faces in zip(bedges, faces_adj_by_edges[bool_boundary_edges]):
        centroids_of_nodes_edge = nodes_centroids[nodes_of_edges[edge]]
        face_to_edge_distance[edge][0] = distance_from_point_to_line(faces_centroids[faces[0]], centroids_of_nodes_edge[0], centroids_of_nodes_edge[1])
    
    return face_to_edge_distance
    
def ordenate_nodes_of_edges(edges, faces_adj_by_edges, nodes_of_faces, nodes_of_edges):
    
    for edge, faces in zip(edges, faces_adj_by_edges):
        nodes_of_face0 = nodes_of_faces[faces[0]]
        nedges = nodes_of_edges[edge]
        position_1 = np.argwhere(nodes_of_face0==nedges[0]).flatten()[0]
        position_2 = np.argwhere(nodes_of_face0==nedges[1]).flatten()[0]
        n_nodes_of_face0 = len(nodes_of_face0)
        
        if position_1 == n_nodes_of_face0-1 and position_2 == 0:
            pass
        elif position_2 == n_nodes_of_face0-1 and position_1 == 0:
            nodes_of_edges[edge][:] = np.flip(nedges)
        elif position_1 > position_2:
            nodes_of_edges[edge][:] = np.flip(nedges)
        
def ordenate_edges_and_nodes_of_nodes_xy_plane(nodes, edges, nodes_adj_by_nodes, edges_adj_by_nodes, nodes_centroids):
    
    for node in nodes:
        nodes_adj = nodes_adj_by_nodes[node]
        edges_adj = edges_adj_by_nodes[node]
        
        centroids_nodes_adj = nodes_centroids[nodes_adj]
        index_sorted = sort_radial_sweep(centroids_nodes_adj, np.arange(len(centroids_nodes_adj)))
        
        nodes_adj[:] = nodes_adj[index_sorted]
        edges_adj[:] = edges_adj[index_sorted]
        centroids_nodes_adj[:] = nodes_centroids[nodes_adj]
        
        new_index = sort_vertices_by_zdirection_xy_plane(centroids_nodes_adj)
        
        nodes_adj[:] = nodes_adj[new_index]
        edges_adj[:] = edges_adj[new_index]       
        
        nodes_adj_by_nodes[node][:] = nodes_adj
        edges_adj_by_nodes[node][:] = edges_adj
        
def ordenate_faces_of_nodes_xy_plane(faces_centroids, faces_adj_by_nodes, nodes_centroids):
    
    for i, faces in enumerate(faces_adj_by_nodes):
        centroid_faces = faces_centroids[faces][:]
        
        if len(faces) == 2:
            centroid_node = nodes_centroids[i]
            ordenate_index = sort_radial_sweep_by_centroid(centroid_node, centroid_faces, np.arange(len(centroid_faces)))
            centroid_faces[:] = centroid_faces[ordenate_index]
            faces2 = faces[ordenate_index]
            new_index = sort_vertices_by_zdirection_xy_plane_by_centroid(centroid_node, centroid_faces)
        
        elif len(faces) > 2:    
            ordenate_index = sort_radial_sweep(centroid_faces, np.arange(len(centroid_faces)))
            centroid_faces[:] = centroid_faces[ordenate_index]
            faces2 = faces[ordenate_index]
            new_index = sort_vertices_by_zdirection_xy_plane(centroid_faces)
        
        else:
            continue
        
        centroid_faces[:] = centroid_faces[new_index]
        faces_adj_by_nodes[i][:] = faces2[new_index]

def define_bool_boundary_nodes(bool_boundary_edges, nodes_of_edges, nodes):
    
    bool_boundary_nodes = np.full(len(nodes), False, dtype=bool)
    nodes_in_boundary = np.unique(nodes_of_edges[bool_boundary_edges].flatten())
    bool_boundary_nodes[nodes_in_boundary] = True
    return bool_boundary_nodes
    