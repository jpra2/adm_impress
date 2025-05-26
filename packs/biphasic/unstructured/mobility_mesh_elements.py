from packs.manager import BoundaryConditions
from packs.biphasic.mobility import BiphasicMobility
from packs.biphasic.relative_perm.brooks_and_corey import BrooksAndCorey

import numpy as np
from typing import Tuple
import time


def calculate_mean(volumes: np.ndarray, mobility: np.ndarray) -> float:
    # import pdb; pdb.set_trace()
    return (volumes.dot(mobility))/volumes.sum()

def nodes_mobility(total_mobility: np.ndarray, volumes: np.ndarray, nodes_elements_adjacencies: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    
    nnodes = len(nodes_elements_adjacencies)
    mob = np.zeros(nnodes)

    for node in range(nnodes):
        elements = nodes_elements_adjacencies[node]
        areas = volumes[elements]
        mob_elements = total_mobility[elements]
        mob[node] = calculate_mean(areas, mob_elements)
    
    return mob
    
def edges_mobility(total_mobility_nodes: np.ndarray, nodes_of_edges: np.ndarray, bc: BoundaryConditions) -> Tuple[np.ndarray, np.ndarray]:
    mob_nodes_adjs = total_mobility_nodes[nodes_of_edges]

    mob_edges = np.mean(mob_nodes_adjs, axis=1)
    
    return mob_edges

def direct_edges_mobility(total_mobility: np.ndarray, areas: np.ndarray, nodes_elements_adjacencies: np.ndarray, nodes_of_edges: np.ndarray, bc: BoundaryConditions, biphasic_mobility: BiphasicMobility, relative_perm: BrooksAndCorey) -> Tuple[np.ndarray, np.ndarray]:
    mob_nodes = nodes_mobility(total_mobility, areas, nodes_elements_adjacencies)
    mob_edges = edges_mobility(mob_nodes, nodes_of_edges, bc)

    edges_presc_sat = bc['water_saturation_edges']['id']
    if len(edges_presc_sat) > 0:
        saturation_edges = bc['water_saturation_edges']['value']
        krw_edges, kro_edges = relative_perm.calculate(saturation_edges)
        mobw_edges, mobo_edges = biphasic_mobility.calculate(krw_edges, kro_edges)
        total_mobility_edges = biphasic_mobility.get_total_mobility(mobw_edges, mobo_edges)

        mob_edges[edges_presc_sat] = total_mobility_edges
        

    return mob_edges






