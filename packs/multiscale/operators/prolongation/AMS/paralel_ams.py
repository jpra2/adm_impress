import numpy as np
import scipy.sparse as sp
from scipy.sparse import linalg
from .ams_tpfa import AmsTpfa
from .ams_mpfa import AMSMpfa
from .....common_files.common_infos import CommonInfos
import multiprocessing as mp

class InfoForLocalOperator(CommonInfos):

	def __init__(self,
		T: 'transmissibility matrix',
		gids: 'global ids',
		dual_id: 'global dual ids',
		primal_id: 'global primal ids'):

		self.T = T
		self.gids = gids
		self.dual_id = dual_id
		self.primal_id = primal_id


class SubDomain:
	def __init__(self, dual_info: 'gids of dual volume'):
		self.dual_info = dual_info

class LocalOperator:
	def __init__(self, subDomains: 'list of SubDomain',
		tpfa=True):

		self.subdomains = subDomains
		self.info = info
		self.tpfa = tpfa

	def get_local_t(self, T, volumes):
		T2 =  T[volumes][:,volumes]
		data = np.array(T2.sum(axis=1).transpose())[0]
		data2 = T2.diagonal()
		data2 -= data
		T2.setdiag(data2)
		return T2

	def local_operator(self):
		g_dual_id = self.info.dual_id
		g_primal_id = self.info.primal_id
		T = self.info.T
		gids = self.info.gids
		n_gids = len(gids)

		for volumes in self.dual_info:
			interns = volumes[g_dual_id[volumes]==0]
			faces = volumes[g_dual_id[volumes]==1]
			edges = volumes[g_dual_id[volumes]==2]
			vertexes = volumes[g_dual_id[volumes]==3]
			n_volumes = len(volumes)
			local_ids = np.arange(n_volumes)
			primal_ids_vertices = g_primal_id[vertexes]
			primal_ids = g_primal_id[volumes]
			local_primal_ids_vertices = np.arange(len(vertexes))
			map_g_to_local_primal_id = dict(zip(primal_ids_vertices, local_primal_ids))
			map_local_to_g_primal_id = dict(zip(local_primal_ids, primal_ids_vertices))
			local_primal_ids = np.array([map_g_to_local_primal_id[k] for k in primal_ids])
			remap_gids = gids.copy()
			remap_gids[volumes] = local_ids
			interns_local = remap_gids[interns]
			edges_local = remap_gids[edges]
			faces_local = remap_gids[faces]
			vertexes_local = remap_gids[vertexes]
			T2 =  self.get_local_t(T, volumes)

			if self.tpfa:
				operator = AmsTpfa
			else:
				operator = AMSMpfa

			operator = operator(interns_local, edges_local, faces_local, vertexes_local, local_ids, local_primal_ids)
			











		s_gids = set(gids)
		n1 = len(gids)
		local_ids = np.arange(n1)
		primal_ids_local = np.arange(nv)

		primal_ids_vertices = primal_ids_global[vertices]
		remap_primal_ids = np.arange(primal_ids_vertices.max()+1)
		remap_primal_ids[primal_ids_vertices] = primal_ids_local

		nni = ni
		nnf = nf + nni
		nne = nnf + ne
		nnv = nne + nv

		gids2 = gids_global
		remap_gids2 = gids2.copy()
		remap_gids2[gids] = local_ids

		T2 = T[gids][:,gids]
		data = np.array(T2.sum(axis=1).transpose())[0]
		data2 = T2.diagonal()
		data2 -= data
		T2.setdiag(data2)

		primal_ids_local_local_ids = remap_primal_ids[primal_ids_global[gids]]

		amstpfa = AmsTpfa(
			local_ids[0:nni],
			local_ids[nni:nnf],
			local_ids[nnf:nne],
			local_ids[nne:nnv],
			local_ids,
			primal_ids_local_local_ids)

		op_local = amstpfa.run(T2)

		d0 = sp.find(op_local)

    	lines_global_op = gids[d0[0]]
    	cols_global_op = primal_ids[d0[1]]
    	data_global_op = d0[2]

    	return lines_global_op, cols_global_op, data_global_op

if __name__=="__main__":
	nworker=6
	#balanceamento para obter entradas
	subdomains=[SubDomain(entradas[i]) for i in range(nworker)]

	procs = [mp.Process(target=LocalOperator, args=[s]) for s in subDomains] #Cria processos

	for p in procs:
        p.start()

	## receber resultados
	for p in procs:
        p.join()
