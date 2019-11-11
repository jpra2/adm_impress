from.data import Data
from ..directories import variables_loaded

def initDataClass(M):
    n_nodes = len(M.nodes)
    n_faces = len(M.faces)
    n_edges = len(M.edges)
    n_volumes = len(M.volumes)

    data = Data(n_nodes, n_faces, n_edges, n_volumes, M)
    verif_format = ['float', 'int']

    data_variables = variables_loaded

    variables = data_variables.keys()

    for variable in variables:
        try:
            infos = data_variables[variable]
        except:
            continue

        data.get_info_data(variable, infos['data size'], infos['data type'], infos['type'])
