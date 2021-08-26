from packs.data_class.compositional_cumulative_datamanager import CumulativeCompositionalDataManager
import numpy as np
import matplotlib.pyplot as plt

# description = 'case1_finescale_'
description = 'case3_finescale_3k'
# description = 'case2_adm_'
# compositional_data = CompositionalData(description=description)
case1 = CumulativeCompositionalDataManager(description=description)
datas_case1 = case1.load_all_datas_from_keys(['loop_array'])

# description = 'case2_adm_'
# description = 'case4_adm_3k'
# description = 'case5_adm_3k'
description = 'case7_adm_3k'
case2 = CumulativeCompositionalDataManager(description=description)
# datas_case2 = case2.load_all_datas()
datas_case2 = case2.load_all_datas_from_keys(['loop_array'])

n_cases = min(len(datas_case1), len(datas_case2))

loops1 = []
for data in datas_case1:
    loops1.append(data['loop_array']['loop'][0])

loops1 = np.array(loops1)

loops2 = []
for data in datas_case2:
    loops2.append(data['loop_array']['loop'][0])
loops2 = np.array(loops2)

n_cases1 = len(datas_case1)
n_cases2 = len(datas_case2)

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

def get_data_from_loop_array(keyword, data_case):
    data = []
    for sim in data_case:
        data.append(sim['loop_array'][keyword][0])
    
    return np.array(data)

datas_case1 = organize_cases_by_loop(n_cases1, datas_case1, loops1)
datas_case2 = organize_cases_by_loop(n_cases2, datas_case2, loops2)

loops1 = np.sort(loops1)

case1_oil_production = get_data_from_loop_array('oil_production', datas_case1)
case2_oil_production = get_data_from_loop_array('oil_production', datas_case2)

case1_gas_production = get_data_from_loop_array('gas_production', datas_case1)
case2_gas_production = get_data_from_loop_array('gas_production', datas_case2)

case1_time = get_data_from_loop_array('t', datas_case1)
case2_time = get_data_from_loop_array('t', datas_case2)
case1_time = case1_time/86400
case2_time = case2_time/86400

n_volumes_update = get_data_from_loop_array('n_volumes_update_base_functions', datas_case2)
total_volumes_updated = get_data_from_loop_array('total_volumes_updated', datas_case2)
max_time_case2 = case2_time.max()
max_time_case1 = case1_time.max()
max_time = min([case1_time.max(), case2_time.max()])

test1 = case1_time <= max_time
test2 = case2_time <= max_time

case1_time = case1_time[test1]
case1_gas_production = case1_gas_production[test1]
case1_oil_production = case1_oil_production[test1]

case2_time = case2_time[test2]
case2_gas_production = case2_gas_production[test2]
case2_oil_production = case2_oil_production[test2]

n_volumes_update = n_volumes_update[test2]
total_volumes_updated = total_volumes_updated[test2]

# import pdb; pdb.set_trace()


# import pdb; pdb.set_trace()


# p1 = extrair_dado(datas_case1, 'pressure')
# p2 = extrair_dado(datas_case2, 'pressure')

# erro = erro_abs(p1, p2)
# erro_rel = erro / p1

# erro_l2 = np.linalg.norm(erro_rel, axis=1)
# erro_max = np.max(erro_rel, axis=1)


fig, (ax1, ax2) = plt.subplots(2, 1)

ax1.plot(case1_time, case1_oil_production, '-', label='Finescale')
ax1.plot(case2_time, case2_oil_production, '-', label='Adm')
ax1.set_ylabel('Oil production')
ax1.set_xlabel('time [days]')
ax1.legend()

ax2.plot(case1_time, case1_gas_production, '-', label='Finescale')
ax2.plot(case2_time, case2_gas_production, '-', label='Adm')
ax2.set_ylabel('Gas production')
ax2.set_xlabel('time [days]')
# ax2.set_ylabel('Norma $L_{\infty}$')
# ax2.set_xlabel('loops')
ax2.legend()

# plt.subplots_adjust(left=0.1,
#                     bottom=0.1, 
#                     right=0.9, 
#                     top=0.9, 
#                     wspace=0.4, 
#                     hspace=0.4)

fig.tight_layout()

# plt.savefig('figura1.png')

# # fig, axes = plt.subplots(2,2)

# # x_axis = np.arange(len(p1[0]))


# # axes[0,0].plot(x_axis, p1[0], '.-')
# # axes[0,0].set_ylabel('loop 1')
# # axes[0,1].plot(x_axis, p1[50], '.-')
# # axes[0,1].set_ylabel('loop 50')
# # axes[1,0].plot(x_axis, p1[100], '.-')
# # axes[1,0].set_ylabel('loop 100')
# # axes[1,1].plot(x_axis, p1[149], '.-')
# # axes[1,1].set_ylabel('loop 150')

# # fig.suptitle('Campo de pressão')
# # fig.tight_layout()

plt.savefig('figura7.png')

fig, (ax1, ax2) = plt.subplots(2, 1)
ax1.plot(case2_time, n_volumes_update, '-')
ax1.set_xlabel('time [days]')
ax1.set_ylabel('N volumes')

ax2.plot(case2_time, total_volumes_updated, '-')
ax2.set_xlabel('time [days]')
ax2.set_ylabel('Total volumes')


plt.subplots_adjust(left=0.1,
                    bottom=0.1, 
                    right=0.9, 
                    top=0.9, 
                    wspace=0.4, 
                    hspace=0.4)
# fig.tight_layout()
fig.suptitle('Volumes para atualizar as funcoes de base')

plt.savefig('figura8.png')



import pdb; pdb.set_trace()


        
    