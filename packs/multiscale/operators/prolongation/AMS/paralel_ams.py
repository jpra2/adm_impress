import numpy as np
import scipy.sparse as sp
from scipy.sparse import linalg
from .ams_tpfa import AmsTpfa
import multiprocessing as mp

class SubDomain:
	def __init__(self):
		self.interns,
		self.faces,
		self.edges,
		self.vertices,
		self.gids_global,
		self.primal_ids_global,
		self.T: 'transmissibility matrix'

class LocalOperator
	def __init__(subDomain):
		self.local_OP=self.local_operator(subDomain)

	def local_operator(subDomain):
		interns=subDomain.intens

		gids = np.concatenate([interns, faces, edges, vertices])
		ni = len(interns)
		nf = len(faces)
		ne = len(edges)
		nv = len(vertices)
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
