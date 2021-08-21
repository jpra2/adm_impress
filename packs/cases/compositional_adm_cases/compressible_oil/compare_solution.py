from packs.data_class.compositional_cumulative_datamanager import CumulativeCompositionalDataManager
import numpy as np
import matplotlib.pyplot as plt

description = 'case1_finescale_'
# description = 'case2_adm_'
# compositional_data = CompositionalData(description=description)
case1 = CumulativeCompositionalDataManager(description=description)
datas_case1 = case1.load_all_datas()

description = 'case2_adm_'
case2 = CumulativeCompositionalDataManager(description=description)
datas_case2 = case2.load_all_datas()

n_cases = min(len(datas_case1), len(datas_case2))

loops1 = []
for data in datas_case1:
    loops1.append(data['loop_array']['loop'][0])

loops1 = np.array(loops1)

loops2 = []
for data in datas_case2:
    loops2.append(data['loop_array']['loop'][0])
loops2 = np.array(loops2)

def organize_cases_by_loop(n_cases, datas, loops):
    resp = []
    for i in range(1, n_cases+1):
        try:
            resp.append(datas[np.argwhere(loops == i)[0][0]])
        except:
            import pdb; pdb.set_trace()
    
    return resp

def extrair_dado(data_case, key):
    resp = []
    for data in data_case:
        resp.append(data[key])
    
    resp = np.array(resp)
    return resp

def erro_abs(v1, v2):
    return np.absolute(v1 - v2)

datas_case1 = organize_cases_by_loop(n_cases, datas_case1, loops1)
datas_case2 = organize_cases_by_loop(n_cases, datas_case2, loops2)

loops1 = np.sort(loops1)

p1 = extrair_dado(datas_case1, 'pressure')
p2 = extrair_dado(datas_case2, 'pressure')

erro = erro_abs(p1, p2)
erro_rel = erro / p1

erro_l2 = np.linalg.norm(erro_rel, axis=1)
erro_max = np.max(erro_rel, axis=1)


# fig, (ax1, ax2) = plt.subplots(2, 1)

# ax1.plot(loops1, erro_l2, '.-')
# ax1.set_ylabel('Norma L2')

# ax2.plot(loops1, erro_max, '.-')
# ax2.set_ylabel('Norma $L_{\infty}$')
# ax2.set_xlabel('loops')

fig, axes = plt.subplots(2,2)

x_axis = np.arange(len(p1[0]))


axes[0,0].plot(x_axis, p1[0], '.-')
axes[0,0].set_ylabel('loop 1')
axes[0,1].plot(x_axis, p1[50], '.-')
axes[0,1].set_ylabel('loop 50')
axes[1,0].plot(x_axis, p1[100], '.-')
axes[1,0].set_ylabel('loop 100')
axes[1,1].plot(x_axis, p1[149], '.-')
axes[1,1].set_ylabel('loop 150')

fig.suptitle('Campo de pressão')
fig.tight_layout()

plt.savefig('figura2.png')





import pdb; pdb.set_trace()


        
    