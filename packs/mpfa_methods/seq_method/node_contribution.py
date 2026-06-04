from packs.biphasic.relative_perm.brooks_and_corey import BrooksAndCorey
from packs.biphasic.mobility import BiphasicMobility
from packs.biphasic.unstructured.mobility_mesh_elements import direct_edges_mobility
from packs.manager import MeshProperty, BoundaryConditions

from packs.examples.same_functions import (
    update_xi_params
)


import numpy as np
from typing import Any



def calculate_total_mobility_by_node(
    relative_perm: BrooksAndCorey,
    biphasic_mobility: BiphasicMobility,
    s: np.ndarray,
    faces_node
) -> np.ndarray:
    krw_faces_node, kro_faces_node = relative_perm.calculate(s[faces_node])
    mobw_faces_node, mobo_faces_node = biphasic_mobility.calculate(krw_faces_node, kro_faces_node)
    total_mobility_faces_node = biphasic_mobility.get_total_mobility(mobw_faces_node, mobo_faces_node)
    return total_mobility_faces_node

def update_total_mobility_edges(
    relative_perm: BrooksAndCorey,
    biphasic_mobility: BiphasicMobility,
    saturation: np.ndarray,
    fp: MeshProperty,
    bc: BoundaryConditions
):
    krw_faces, kro_faces = relative_perm.calculate(saturation)
    mobw_faces, mobo_faces = biphasic_mobility.calculate(krw_faces, kro_faces)
    total_mobility_faces = biphasic_mobility.get_total_mobility(mobw_faces, mobo_faces)
    fw_faces = biphasic_mobility.get_fw(mobw_faces, mobo_faces)
    total_mobility_edges = direct_edges_mobility(
        total_mobility_faces,
        fp['areas'],
        fp['faces_of_nodes'],
        fp['nodes_of_edges'],
        bc,
        biphasic_mobility,
        relative_perm
    )

    fp.insert_or_update_data({
        'edges_multiplier': total_mobility_edges
    })

    fp.insert_or_update_data({
        'xi_params': update_xi_params(fp['xi_params_backup'], total_mobility_edges)
    })

def calculate_node_contribution(
    node: np.ndarray[tuple[int], np.dtype[np.uint64]],
    # edges_flux: np.ndarray[tuple[int], np.dtype[np.float64]],
    nodes_weights: np.ndarray,
    xi_params: np.ndarray,
    edges_of_nodes: np.ndarray,
    nodes_of_edges: np.ndarray,
    adjacencies: np.ndarray
):
    
    
    edges_node = edges_of_nodes[node]
    nodes_edges_node = nodes_of_edges[edges_node]
    faces_node = nodes_weights['face_id'][nodes_weights['node_id']==node]
    weights_node = nodes_weights['weight'][nodes_weights['node_id']==node]
    adjacencies_edges_node = adjacencies[edges_node]
    xi_params_edges_node = xi_params[edges_node]
    
    A_node = nodes_edges_node[:, 1] == node
    B_node = nodes_edges_node[:, 0] == node
    
    K_faces = adjacencies_edges_node[:, 0]
    L_faces = adjacencies_edges_node[:, 1]
    test_l_faces = L_faces != -1
    test_k_faces = K_faces != -1
    xi_A = xi_params_edges_node[:, 2]
    xi_B = xi_params_edges_node[:, 3]
        
        
    
    
    
    
    


def mount_edge_flux_internal(
    edge: int,
    weights: np.ndarray,
    adjacencies: np.ndarray,
    eta_params: np.ndarray,
    nodes_of_edges: np.ndarray,
    p: np.ndarray,
    s: np.ndarray,
    area: np.ndarray,
    biphasic_mobility: BiphasicMobility,
    relative_perm: BrooksAndCorey
):
    Anode, Bnode = nodes_of_edges[edge]
    Kface, Lface = adjacencies[edge]
    
    weights_Anode = weights[weights['node_id'] == Anode]
    weights_Bnode = weights[weights['node_id'] == Bnode]
    
    faces_Anode = weights_Anode['face_id']
    faces_Bnode = weights_Bnode['face_id']
    
    total_mobility_faces_Anode = calculate_total_mobility_by_node(
        relative_perm,
        biphasic_mobility,
        s,
        faces_Anode
    )
    
    total_mobility_faces_Bnode = calculate_total_mobility_by_node(
        relative_perm,
        biphasic_mobility,
        s,
        faces_Bnode
    )
    
    mA = (total_mobility_faces_Anode*area[faces_Anode]).sum()/area[faces_Anode].sum()
    mB = (total_mobility_faces_Bnode*area[faces_Bnode]).sum()/area[faces_Bnode].sum()
    
    total_mobility_edge = (mA + mB)/2
    
    pA = (weights_Anode['weight']*p[weights_Anode['face_id']]).sum()
    pB = (weights_Bnode['weight']*p[weights_Bnode['face_id']]).sum()
    
    edge_flux = total_mobility_edge*(eta_params*[p[Kface], p[Lface], pB, pA]).sum()
    
    return edge_flux

def mount_residual(
    p: np.ndarray,
    s: np.ndarray,
    eta_params: np.ndarray,
    adjacencies: np.ndarray,
    weights: np.ndarray,
    faces_of_nodes: np.ndarray,
    nodes_of_edges: np.ndarray,
    edges_dim: np.ndarray,
    area: np.ndarray,
    edges: np.ndarray,
    bool_internal_edges: np.ndarray
):
    
    internal_edges = edges[bool_internal_edges]
    edges_flux = 0
    
    
    
    
    