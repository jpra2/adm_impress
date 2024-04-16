import numpy as np
import scipy.sparse as sp
from packs import defnames


class Adm:

    def update_levels(
            self,
            list_primal_ids,
            fine_levels,
            fine_ids,
            n_levels,
            fine_faces_of_faces
    ):
        
        fine_ids_prox_level = fine_ids[fine_levels == 0]

        for level in range(1, n_levels):
            primal_ids = list_primal_ids[level]
            test_level_ant = fine_levels == level - 1
            fine_coarse_ids_ant = np.unique(primal_ids[test_level_ant])
            level_ant_bool = np.isin(primal_ids, fine_coarse_ids_ant)
            fine_levels[level_ant_bool] = level - 1
            fine_ids_level0 = fine_ids[fine_levels <= level - 1]
            fine_ids_prox_level = np.setdiff1d(
                np.unique(np.concatenate(fine_faces_of_faces[fine_ids_level0])),
                fine_ids_level0
            )
            fine_levels[fine_ids_prox_level] = level
        
        not_modify = (fine_levels <= level - 1) & (fine_levels != -1)
        modify = ~not_modify
        fine_levels[modify] = level

    def get_adm_prolongation_operator(
            self,
            list_OP,
            list_OR,
            fine_levels: np.ndarray,
            list_primal_ids,
            list_dual_ids,
            fine_ids,
            n_levels
    ):
        if n_levels > 2:
            raise NotImplementedError
        
        list_ids_level = []
        for level in range(n_levels):
            list_ids_level.append(np.unique(list_primal_ids[level][fine_levels==level]))
        
        ### do nivel 0 para o nivel 1
        level = 1
        all_n_adm_ids = [list_primal_ids[0].shape[0]]
        all_adm_ids = [list_primal_ids[0]]

        dual_ids_level_ant = list_dual_ids[level-1]
        
        OP = list_OP[level-1]
        OR = list_OR[level-1]
        data_op = sp.find(OP)
        data_or = sp.find(OR)

        primal_ids_ant = list_primal_ids[level-1]
        primal_ids_prox = list_primal_ids[level]
        fine_vertices = primal_ids_ant[dual_ids_level_ant == defnames.dual_ids('vertice_id')]
        coarse_ids_fine_vertices = primal_ids_prox[fine_vertices]

        ids_level_prox = np.unique(primal_ids_prox[fine_levels >= level])
        coarse_ids_to_remove = np.unique(primal_ids_prox[fine_levels < level])
        ids_level_ant = np.unique(primal_ids_ant[fine_levels < level])

        n_fine = ids_level_ant.shape[0] ## os que se mantem no nivel
        n_coarse = ids_level_prox.shape[0] ## os que foram engrossados

        n_coarse_adm = n_fine + n_coarse
        adm_ids = np.arange(n_coarse_adm)
        all_n_adm_ids.append(adm_ids.shape[0])
        all_adm_ids.append(adm_ids)

        ## removendo dos operadores as linhas e as colunas que vao ficar no nivel anterior
        test1 = np.isin(data_op[0], ids_level_ant)
        test = ~test1
        lines_op_adm = data_op[0][test]
        cols_op_adm = data_op[1][test]
        data_op_adm = data_op[2][test]

        test5 = np.isin(data_or[1], ids_level_ant)
        test5 = ~test5
        lines_or_adm = data_or[0][test5]
        cols_or_adm = data_or[1][test5]
        data_or_adm = data_or[2][test5]

        test2 = np.isin(cols_op_adm, coarse_ids_to_remove)
        coarse_ids_to_modify_by_adm_vertice_id = np.unique(cols_op_adm[test2])
        test6 = np.isin(lines_or_adm, coarse_ids_to_remove)
        or_coarse_ids_to_modify_by_vertice = np.unique(lines_or_adm[test6])

        ### criando o remapeamento adm
        ## os que permanecem no nivel
        remap_fine_cols = np.repeat(-1, primal_ids_ant.shape[0])
        remap_fine_cols[ids_level_ant] = np.arange(n_fine)
        ## os que vao para o proximo nivel
        remap_coarse_cols = np.repeat(-1, primal_ids_prox.shape[0])
        remap_coarse_cols[ids_level_prox] = np.arange(n_fine, n_fine + n_coarse)

        ## atualizando os valores do remapeamento adm
        test3 = ~test2
        cols_op_adm[test3] = remap_coarse_cols[cols_op_adm[test3]]
        for cid in coarse_ids_to_modify_by_adm_vertice_id:
            fine_vertice = fine_vertices[coarse_ids_fine_vertices==cid]
            test4 = cols_op_adm == cid
            cols_op_adm[test4] = remap_fine_cols[fine_vertice]
        
        test7 = ~test6
        lines_or_adm[test7] = remap_coarse_cols[lines_or_adm[test7]]
        for cid in or_coarse_ids_to_modify_by_vertice:
            fine_vertice = fine_vertices[coarse_ids_fine_vertices==cid]
            test8 = lines_or_adm == cid
            lines_or_adm[test8] = remap_fine_cols[fine_vertice]

        ### adicionando os que permanecem no nivel
        lines_op_adm = np.append(lines_op_adm, ids_level_ant)
        cols_op_adm = np.append(cols_op_adm, remap_fine_cols[ids_level_ant])
        data_op_adm = np.append(data_op_adm, np.repeat(1.0, ids_level_ant.shape[0]))

        lines_or_adm = np.append(lines_or_adm, remap_fine_cols[ids_level_ant])
        cols_or_adm = np.append(cols_or_adm, ids_level_ant)
        data_or_adm = np.append(data_or_adm, np.repeat(1.0, ids_level_ant.shape[0]))

        OP_adm = sp.csc_matrix((data_op_adm, (lines_op_adm, cols_op_adm)), shape=(all_n_adm_ids[level-1], all_n_adm_ids[level]))
        OR_adm = sp.csc_matrix((data_or_adm, (lines_or_adm, cols_or_adm)), shape=(all_n_adm_ids[level], all_n_adm_ids[level-1]))

        return OP_adm, OR_adm

