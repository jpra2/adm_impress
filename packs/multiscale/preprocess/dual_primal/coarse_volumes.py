

class CoarseVolume:

    def __init__(self, volumes=None, faces=None, edges=None, nodes=None,
                 boundary_faces=None, primal_id=None, element=None, volumes_element=None,
                 faces_element=None, edges_element=None, nodes_element=None):

        self.volumes = volumes
        self.faces = faces
        self.edges = edges
        self.nodes = nodes
        self.boundary_faces = boundary_faces
        self.primal_id = primal_id
        self.element = element

        self.volumes_element = volumes_element
        self.faces_element = faces_element
        self.edges_element = edges_element
        self.nodes_element = nodes_element

