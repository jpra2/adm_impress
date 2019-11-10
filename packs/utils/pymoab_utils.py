import numpy as np
from pymoab import core, types, rng, topo_util, skinner
import os
import yaml
import sys


def get_faces(mb, elements):
    mtu = topo_util.MeshTopoUtil(mb)
    elements = get_elements(mb, elements)
    faces = mtu.get_bridge_adjacencies(elements, 3, 2)
    return faces

def get_boundary_of_volumes(mb, elements):
    """
    input:
        mb: core of pymoab
        elements: meshset or elements of the mesh
    output:
        boundary of meshset or elements
    """
    faces = UtilsPymoab.get_faces(mb, elements)

    bound_faces = []

    for face in faces:
        elems = mb.get_adjacencies(face, 3)
        if len(elems) < 2:
            bound_faces.append(face)
        elif elems[0] in elements and elems[1] not in elements:
            bound_faces.append(face)
        elif elems[1] in elements and elems[0] not in elements:
            bound_faces.append(face)

    return rng.Range(bound_faces)

def get_faces_in_intersection_between_volumes(mb, elements1, elements2):
    """
    Retorna as faces na interseccao entre dois ranges de volumes ou dois meshsets
    """
    bound_faces1 = UtilsPymoab.get_boundary_of_volumes(mb, elements1)
    bound_faces2 = UtilsPymoab.get_boundary_of_volumes(mb, elements2)
    return rng.intersect(bound_faces1, bound_faces2)

def get_elements(mb, elements):
    """
    retorna o rng.Range dos elementos de entrada
    """
    if isinstance(elements, rng.Range):
        return elements
    elif isinstance(elements, int):
        elements = mb.get_entities_by_handle(elements)
        return elements
    elif isinstance(elements, np.ndarray) or isinstance(elements, list):
        return rng.Range(elements)
    else:
        raise ValueError('tipo de dado incorreto')

def get_all_entities(mb):
    mtu = topo_util.MeshTopoUtil(mb)
    root_set = mb.get_root_set()
    all_volumes = mb.get_entities_by_dimension(0, 3)
    all_nodes = mb.get_entities_by_dimension(0, 0)
    mtu.construct_aentities(all_nodes)
    all_faces = mb.get_entities_by_dimension(0, 2)
    all_edges = mb.get_entities_by_dimension(0, 1)

    return [all_nodes, all_edges, all_faces, all_volumes]

def get_all_tags_1(mb, list_names_tags):
    tags = []
    for name in list_names_tags:
        tag = mb.tag_get_handle(name)
        tags.append(tag)

    return tags

def get_all_tags_2(mb, list_names_tags):
    tags = {}
    for name in list_names_tags:
        try:
            tag = mb.tag_get_handle(str(name))
        except:

            print(name, 'Nao existe no arquivo')
            continue
            # sys.exit(0)
            # import pdb; pdb.set_trace()

        tags[name] = tag

    return tags

def set_faces_in_boundary_by_meshsets_dep0(mb, mtu, meshsets, faces_boundary_meshset_tag):
    all_faces_on_boundary = mb.create_meshset()
    meshsets = list(meshsets)
    n = len(meshsets)

    ms0 = set()

    # for m1 in meshsets:
    #     ms0.add(m1)
    #     cont = 0
    #     elems1 = mb.get_entities_by_handle(m1)
    #     faces1 = mtu.get_bridge_adjacencies(elems1, 2, 2)
    #     for m2 in meshsets:
    #         if m2 in ms0:
    #             continue
    #         elems2 = mb.get_entities_by_handle(m2)
    #         faces2 = mtu.get_bridge_adjacencies(elems2, 2, 2)
    #         intersect = rng.intersect(faces1, faces2)
    #         if len(intersect) < 1:
    #             continue
    #         cont+=1
    #         mb.add_entities(all_faces_on_boundary, intersect)
    #     print(cont)
    #     if cont > 3:
    #         print(i)
    #         print(cont)
    #         import pdb; pdb.set_trace()
    # import pdb; pdb.set_trace()

    for m1 in meshsets:
        cont = 0
        elems1 = mb.get_entities_by_handle(m1)
        faces1 = mtu.get_bridge_adjacencies(elems1, 3, 2)
        # for face in faces1:
        #     elems = mb.get_adjacencies(face, 3)
        #     if len(elems) < 2:
        #         continue
        #     if elems[0] in elems1 and elems[1] in elems1:
        #         continue
        #     mb.add_entities(all_faces_on_boundary, rng.Range(face))
        elems2 = mtu.get_bridge_adjacencies(elems1, 2, 3)
        elems3 = rng.subtract(elems2, elems1)
        faces3 = mtu.get_bridge_adjacencies(elems3, 3, 2)
        faces_cont = rng.intersect(faces3, faces1)

        mb.add_entities(all_faces_on_boundary, faces_cont)

    mb.tag_set_data(faces_boundary_meshset_tag, 0, all_faces_on_boundary)

def set_faces_in_boundary_by_meshsets(mb, mtu, meshsets, faces_boundary_meshset_tag, M):
    all_faces_on_boundary = mb.create_meshset()
    meshsets = list(meshsets)
    n = len(meshsets)

    ms0 = set()

    # for m1 in meshsets:
    #     ms0.add(m1)
    #     cont = 0
    #     elems1 = mb.get_entities_by_handle(m1)
    #     faces1 = mtu.get_bridge_adjacencies(elems1, 2, 2)
    #     for m2 in meshsets:
    #         if m2 in ms0:
    #             continue
    #         elems2 = mb.get_entities_by_handle(m2)
    #         faces2 = mtu.get_bridge_adjacencies(elems2, 2, 2)
    #         intersect = rng.intersect(faces1, faces2)
    #         if len(intersect) < 1:
    #             continue
    #         cont+=1
    #         mb.add_entities(all_faces_on_boundary, intersect)
    #     print(cont)
    #     if cont > 3:
    #         print(i)
    #         print(cont)
    #         import pdb; pdb.set_trace()
    # import pdb; pdb.set_trace()

    for m1 in meshsets:
        cont = 0
        elems1 = mb.get_entities_by_handle(m1)
        faces1 = mtu.get_bridge_adjacencies(elems1, 3, 2)
        # for face in faces1:
        #     elems = mb.get_adjacencies(face, 3)
        #     if len(elems) < 2:
        #         continue
        #     if elems[0] in elems1 and elems[1] in elems1:
        #         continue
        #     mb.add_entities(all_faces_on_boundary, rng.Range(face))
        elems2 = mtu.get_bridge_adjacencies(elems1, 2, 3)
        elems3 = rng.subtract(elems2, elems1)
        faces3 = mtu.get_bridge_adjacencies(elems3, 3, 2)
        faces_cont = rng.intersect(faces3, faces1)

        mb.add_entities(all_faces_on_boundary, faces_cont)

    mb.tag_set_data(faces_boundary_meshset_tag, M.core.root_set, all_faces_on_boundary)

def load_adm_mesh():
    n_levels = 3
    parent_dir = os.path.dirname(os.path.abspath(__file__))
    parent_parent_dir = os.path.dirname(parent_dir)
    input_dir = os.path.join(parent_parent_dir, 'input')
    flying_dir = os.path.join(parent_parent_dir, 'flying')
    bifasico_dir = os.path.join(flying_dir, 'bifasico')
    utils_dir = os.path.join(parent_parent_dir, 'utils')
    processor_dir = os.path.join(parent_parent_dir, 'processor')

    os.chdir(input_dir)
    with open("inputs.yaml", 'r') as stream:
        data_loaded = yaml.load(stream)
        # data_loaded = yaml.load(stream, Loader=yaml.FullLoader)
        # data_loaded = yaml.full_load(stream)

    input_file = data_loaded['input_file']
    ext_h5m_adm = input_file + '_malha_adm.h5m'
    ADM = data_loaded['ADM']
    tempos_impr = data_loaded['tempos_vpi_impressao']
    contar_loop = data_loaded['contar_loop']
    contar_tempo = data_loaded['contar_tempo']
    imprimir_sempre = data_loaded['imprimir_sempre']

    mb = core.Core()
    mtu = topo_util.MeshTopoUtil(mb)
    os.chdir(flying_dir)
    mb.load_file(ext_h5m_adm)
    list_names_tags = np.load('list_names_tags.npy')
    # list_names_tags = np.delete(list_names_tags, np.where(list_names_tags == 'L_TOT')[0])
    # names_tags_with_level = np.load('names_tags_with_level.npy')
    tags_1 = get_all_tags_2(mb, list_names_tags)
    os.chdir(parent_dir)

    return mb, mtu, tags_1, input_file, ADM, tempos_impr, contar_loop, contar_tempo, imprimir_sempre

def enumerar_volumes_nivel(mb, meshsets, level):
    name_tag = 'IDS_NA_PRIMAL_' + str(level)
    ids_na_primal_tag = mb.tag_get_handle(name_tag, 1, types.MB_TYPE_INTEGER, types.MB_TAG_SPARSE, True)
    conts = 0
    for m in meshsets:
        elems = mb.get_entities_by_handle(m)
        n = len(elems)
        num1 = conts
        num2 = num1 + n
        mb.tag_set_data(ids_na_primal_tag, elems, np.arange(num1, num2))
        conts += n

    return name_tag, ids_na_primal_tag
