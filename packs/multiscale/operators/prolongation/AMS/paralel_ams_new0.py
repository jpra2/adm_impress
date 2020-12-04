import numpy as np
import scipy.sparse as sp
from scipy.sparse import linalg
from .ams_tpfa_new0 import AMSTpfa
from .....common_files.common_infos import CommonInfos
import multiprocessing as mp
import collections
import time

def run_thing(local_operator_obj):
	local_operator_obj.run()


class SubDomain(CommonInfos):

	def __init__(self,
		T: 'global transmissibility matrix',
		g_local_ids: 'global id of local volumes',
		dual_flags: 'dual flags of g_local_ids',
		primal_ids: 'primal ids of g_local_ids',
		get_correction_term=False,
		B_matrix=None,
		Eps_matrix=None,
		total_source_term=None):

		dual_ids = dual_flags
		self.get_correction_term = get_correction_term

		if get_correction_term:
			self.l_B_matrix = self.get_local_matrix(B_matrix, g_local_ids)
			self.l_Eps_matrix = self.get_local_matrix(Eps_matrix, g_local_ids)
			self.l_total_source_term = total_source_term[g_local_ids]
		else:
			self.l_B_matrix = None
			self.l_Eps_matrix = None
			self.l_total_source_term = None

		self.l_dual_flags = dual_flags
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
		# map_gids = dict(zip(g_local_ids, self.local_ids))
		# map_gids = np.repeat(-1, g_local_ids.max()+1)
		# map_gids[g_local_ids] = self.local_ids
		self.l_primal_ids = np.array([map_primal_ids[k] for k in global_primal_ids])
		# self.l_interns = self.local_ids[local_dual_id==0]
		# self.l_faces = self.local_ids[local_dual_id==1]
		# self.l_edges = self.local_ids[local_dual_id==2]
		# self.l_vertices = self.local_ids[local_dual_id==3]


class LocalOperator:
	def __init__(self,
		subDomains: 'list of SubDomains',
		comm: 'comunicator'
	):

		self.subdomains = subDomains
		self.comm = comm

	def run(self):

		lines = []
		cols = []
		data = []
		all_pcorr = []

		for subd in self.subdomains:
			operator = AMSTpfa(subd.local_ids, subd.l_dual_flags, subd.l_primal_ids, get_correction_term=subd.get_correction_term)

			op, pcorr = operator.run(subd.l_T, total_source_term=subd.l_total_source_term, B_matrix=subd.l_B_matrix, Eps_matrix=subd.l_Eps_matrix)
			subd.local_ids[:] = subd.g_local_ids
			ff = sp.find(op)

			ls = subd.local_ids[ff[0]]
			cs = np.array([subd.r_map_primal_ids[k] for k in ff[1]])
			ds = ff[2]

			lines.append(ls)
			cols.append(cs)
			data.append(ds)
			all_pcorr.append(np.array([subd.g_local_ids, pcorr]))

		structure = np.array([lines, cols, data, all_pcorr])

		self.comm.send(structure)


class MasterOP:
	def __init__(self,
		T: 'global transmissibility matrix',
		all_dual_subsets: 'all dual volumes',
		level,
		get_correction_term=False,
		total_source_term=None,
		B_matrix=None,
		Eps_matrix=None):

		n_cpu = mp.cpu_count()
		self.n_workers = n_cpu
		self.level = level
		self.n_total = T.shape[0]
		self.get_correction_term = get_correction_term

		list_of_subdomains = self.get_list_subdomains(T, all_dual_subsets, get_correction_term, B_matrix=B_matrix,
		 	Eps_matrix=Eps_matrix, total_source_term=total_source_term)
		list_of_process_per_cpu = []
		n_subdomains = len(list_of_subdomains)
		resto = n_subdomains % self.n_workers
		n_process_per_cpu = n_subdomains//self.n_workers

		if n_process_per_cpu > 0:

			for i in range(self.n_workers):
				list_of_process_per_cpu.append(list_of_subdomains[i*n_process_per_cpu:n_process_per_cpu*(i+1)])

			if resto != 0:
				for i in range(resto):
					list_of_process_per_cpu[i].append(list_of_subdomains[-i])

		else:
			self.n_workers = resto

			for i in range(resto):
				list_of_process_per_cpu[i].append(list_of_subdomains[-i])

		self.list_of_process_per_cpu = list_of_process_per_cpu

	def get_list_subdomains(self, T, all_dual_subsets, get_correction_term=False, B_matrix=None, Eps_matrix=None, total_source_term=None):
		list_of_subdomains = []
		for dual_subset in all_dual_subsets:
			sarray = dual_subset
			volumes = sarray['volumes']
			dual_ids1 = sarray['dual_id']
			primal_ids1 = sarray['primal_id']

			# all_edges = volumes[dual_ids1==2]
			# contador = collections.Counter(all_edges)
			# coupled_edges = np.array([k for k, v in contador.items() if v > 1])
			# local_gids = np.unique(volumes)
			local_gids = volumes
			# dual_ids = np.concatenate([dual_ids1[volumes==k] for k in local_gids])
			dual_ids = dual_ids1
			# primal_ids = np.concatenate([primal_ids1[volumes==k] for k in local_gids])
			primal_ids = primal_ids1
			list_of_subdomains.append(SubDomain(T, local_gids, dual_ids, primal_ids, get_correction_term=get_correction_term, B_matrix=B_matrix, Eps_matrix=Eps_matrix, total_source_term=total_source_term))

		return list_of_subdomains

	def get_data_for_op(self, all_lines, all_cols, all_datas, set_lines):

		lines = []
		cols = []
		data = []
		ps = []

		for ls, cs, ds in zip(all_lines, all_cols, all_datas):
			resp = set(ls) - set_lines
			if resp:
				conj_lines = []
				for k in resp:
					indice = np.where(ls == k)[0]
					conj_lines.append(indice)

				conj_lines = np.concatenate(conj_lines)

				resp_ls = ls[conj_lines]
				resp_cols = cs[conj_lines]
				resp_data = ds[conj_lines]

				lines.append(resp_ls)
				cols.append(resp_cols)
				data.append(resp_data)

				set_lines = set_lines | resp

		lines = np.concatenate(lines).astype(np.int64)
		cols = np.concatenate(cols).astype(np.int64)
		data = np.concatenate(data)



		return lines, cols, data, set_lines

	def run(self):

		master2worker = [mp.Pipe() for _ in range(self.n_workers)]
		m2w, w2m = list(zip(*master2worker))
		procs = [mp.Process(target=run_thing, args=[LocalOperator(obj, comm)]) for obj, comm in zip(self.list_of_process_per_cpu, w2m)]
		del self.list_of_process_per_cpu

		lines = []
		cols = []
		data = []
		pcorr = np.zeros(self.n_total)

		set_lines = set()
		for proc in procs:
			proc.start()

		for comm in m2w:
			msg = comm.recv()
			resp = msg
			all_lines = resp[0]
			all_cols = resp[1]
			all_datas = resp[2]
			all_pcorr = resp[3]

			ls, cs, ds, set_lines = self.get_data_for_op(all_lines, all_cols, all_datas, set_lines)

			lines.append(ls)
			cols.append(cs)
			data.append(ds)

			for ppp in all_pcorr:
				pcorr[ppp[0].astype(int)] = ppp[1]

		for proc in procs:
			proc.join()

		lines = np.concatenate(lines).astype(np.int64)
		cols = np.concatenate(cols).astype(np.int64)
		data = np.concatenate(data)

		n_volumes = lines.max() + 1
		n_c_volumes = cols.max() + 1

		OP = sp.csc_matrix((data, (lines, cols)), shape=(n_volumes, n_c_volumes))

		return OP, pcorr
