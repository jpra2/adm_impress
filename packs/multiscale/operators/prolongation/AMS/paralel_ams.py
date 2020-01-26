import numpy as np
import scipy.sparse as sp
from scipy.sparse import linalg
from .ams_tpfa import AMSTpfa
from .ams_mpfa import AMSMpfa
from .....common_files.common_infos import CommonInfos
import multiprocessing as mp
import collections

def run_thing(local_operator_obj):
	local_operator_obj.run()

class SubDomain(CommonInfos):

	def __init__(self,
		T: 'global transmissibility matrix',
		g_local_ids: 'global id of local volumes',
		dual_ids: 'dual id of g_local_ids',
		primal_ids: 'primal ids of g_local_ids',
		coupled_edges: 'coupled edges'=[]
	):

		self.l_T = self.get_local_t(T, g_local_ids) # local transmissibility
		self.g_local_ids = g_local_ids
		local_dual_id = dual_ids
		g_vertices = g_local_ids[local_dual_id==3]
		global_primal_ids = primal_ids

		local_primal_ids = np.arange(len(g_vertices))
		primal_ids_vertices = primal_ids[local_dual_id==3]
		map_primal_ids = dict(zip(primal_ids_vertices, local_primal_ids))
		self.r_map_primal_ids = dict(zip(local_primal_ids, primal_ids_vertices))

		self.local_ids = np.arange(len(g_local_ids))
		map_gids = dict(zip(g_local_ids, self.local_ids))
		self.l_primal_ids = np.array([map_primal_ids[k] for k in global_primal_ids])
		self.l_interns = self.local_ids[local_dual_id==0]
		self.l_faces = self.local_ids[local_dual_id==1]
		self.l_edges = self.local_ids[local_dual_id==2]
		self.l_vertices = self.local_ids[local_dual_id==3]
		self.l_coupled_edges = np.array([map_gids[k] for k in coupled_edges])


class LocalOperator:
	def __init__(self,
		subDomains: 'list of SubDomains',
		comm: 'comunicator'
	):

		self.subdomains = subDomains
		self.comm = comm

	def local_operator_dep0(self):
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

	def run(self):

		lines = []
		cols = []
		data = []

		for subd in self.subdomains:
			operator = AMSTpfa(subd.l_interns, subd.l_faces, subd.l_edges, subd.l_vertices,
				subd.local_ids, subd.l_primal_ids, coupled_edges=subd.l_coupled_edges
			)

			op = operator.run(subd.l_T)
			subd.local_ids[:] = subd.g_local_ids
			ff = sp.find(op)

			ls = subd.local_ids[ff[0]]
			cs = np.array([subd.r_map_primal_ids[k] for k in ff[1]])
			ds = ff[2]

			lines.append(ls)
			cols.append(cs)
			data.append(ds)

		structure = np.array([lines, cols, data])

		self.comm.send(structure)

class MasterOP:
	def __init__(self,
		T: 'global transmissibility matrix',
		all_dual_subsets: 'all dual volumes'
	):

		n_cpu = mp.cpu_count()
		self.n_workers = n_cpu

		list_of_subdomains = self.get_list_subdomains(T, all_dual_subsets)
		list_of_process_per_cpu = []
		n_subdomains = len(list_of_subdomains)
		resto = n_subdomains % self.n_workers
		n_process_per_cpu = n_subdomains//self.n_workers

		for i in range(self.n_workers):
			list_of_process_per_cpu.append(list_of_subdomains[i*n_process_per_cpu:n_process_per_cpu*(i+1)])

		if resto != 0:
			for i in range(resto):
				list_of_process_per_cpu[i].append(list_of_subdomains[-i])

		self.list_of_process_per_cpu = list_of_process_per_cpu

	def get_list_subdomains(self, T, all_dual_subsets):
		list_of_subdomains = []
		for dual_subset in all_dual_subsets:
			sarray = np.concatenate([dual_subset])
			volumes = sarray['volumes']
			dual_ids1 = sarray['dual_id']
			primal_ids1 = sarray['primal_id']

			all_edges = volumes[dual_ids1==2]
			contador = collections.Counter(all_edges)
			coupled_edges = np.array([k for k, v in contador.items() if v > 1])
			local_gids = np.unique(volumes)
			dual_ids = np.concatenate([dual_ids1[volumes==k] for k in local_gids])
			primal_ids = np.concatenate([primal_ids1[volumes==k] for k in local_gids])
			list_of_subdomains.append(SubDomain(T, local_gids, dual_ids, primal_ids, coupled_edges))

		return list_of_subdomains

	def run(self):

		master2worker = [mp.Pipe() for _ in range(self.n_workers)]
		m2w, w2m = list(zip(*master2worker))
		procs = [mp.Process(target=run_thing, args=[LocalOperator(obj, comm)]) for obj, comm in zip(self.list_of_process_per_cpu, w2m)]
		del self.list_of_process_per_cpu

		lines = []
		cols = []
		data = []

		set_lines = set()
		for proc in procs:
			proc.start()

		for comm in m2w:
			msg = comm.recv()
			ls = msg[0]
			cs = msg[1]
			ds = msg[2]

			set_ls = set(ls) & set_lines
			set_m_ls = set_ls - set_lines
			positions = []

			for li in set_ls:
				pos = np.where(ls == li)[0]
				positions.append(pos[0])

			ls = np.delete(ls, positions)
			cs = np.delete(cs, positions)
			ds = np.delete(ds, positions)

			lines.append(ls)
			cols.append(cs)
			data.append(ds)

			set_lines = set_lines | set_m_ls

		for proc in procs:
			proc.join()

		lines = np.concatenate(lines)
		cols = np.concatenate(cols)
		data = np.concatenate(data)

		n_volumes = len(np.unique(lines))
		n_c_volumes = len(np.unique(cols))

		OP = sp.csc_matrix((data, (lines, cols)), shape=(n_volumes, n_c_volumes))

		return OP


















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
