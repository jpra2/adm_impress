from packs.manager import BoundaryConditions

import numpy as np
from typing import Tuple


def calculate_mean(volumes: np.ndarray, mobility: np.ndarray) -> float:
    return (volumes.dot(mobility))/volumes.sum()

def nodes_mobility(water_mobility: np.ndarray, oil_mobility: np.ndarray, volumes: np.ndarray, nodes_elements_adjacencies: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    
    nnodes = len(nodes_elements_adjacencies)
    mobw = np.zeros(nnodes)
    mobo = mobw.copy()

    for node in range(nnodes):
        elements = nodes_elements_adjacencies[node]
        areas = volumes[elements]
        mobw_elements = water_mobility[elements]
        mobo_elements = oil_mobility[elements]
        mobw[node] = calculate_mean(areas, mobw_elements)
        mobo[node] = calculate_mean(areas, mobo_elements)
    
    return mobw, mobo
    
def edges_mobility(mobw_nodes: np.ndarray, mobo_nodes: np.ndarray, nodes_of_edges: np.ndarray, bc: BoundaryConditions) -> Tuple[np.ndarray, np.ndarray]:
    mobw_nodes_adjs = mobw_nodes[nodes_of_edges]
    mobo_nodes_adjs = mobo_nodes[nodes_of_edges]

    mobw_edges = np.mean(mobw_nodes_adjs, axis=1)
    mobo_edges = np.mean(mobo_nodes_adjs, axis=1)

    edges_with_sat_presc = bc['water_saturation_edges']['id']
    if len(edges_with_sat_presc) > 0:
        # TODO implementar essa parte
        raise NotImplementedError
    
    return mobw_edges, mobo_edges

def direct_edges_mobility(water_mobility: np.ndarray, oil_mobility: np.ndarray, areas: np.ndarray, nodes_elements_adjacencies: np.ndarray, nodes_of_edges: np.ndarray, bc: BoundaryConditions) -> Tuple[np.ndarray, np.ndarray]:
    mobw_nodes, mobo_nodes = nodes_mobility(water_mobility, oil_mobility, areas, nodes_elements_adjacencies)
    mobw_edges, mobo_edges = edges_mobility(mobw_nodes, mobo_nodes, nodes_of_edges, bc)

    return mobw_edges, mobo_edges






