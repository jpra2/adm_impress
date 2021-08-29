from packs.data_class.compositional_cumulative_datamanager import CumulativeCompositionalDataManager
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# description = 'case1_finescale_'
description1 = 'case3_finescale_3k'
# description = 'case2_adm_'
# compositional_data = CompositionalData(description=description)
case1 = CumulativeCompositionalDataManager(description=description1)
datas_case1 = case1.load_all_datas_from_keys(['loop_array'])

fig_str = 'figura_'

# description = 'case2_adm_'
# description = 'case4_adm_3k'
# description = 'case5_adm_3k'
# description = 'case7_adm_3k'
# description = 'case8_adm_3k'
description2 = 'case9_adm_3k'
# description2 = 'case10_adm_3k'
case2 = CumulativeCompositionalDataManager(description=description2)
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

def organize_cases_by_loop(datas, loops):
    loops2 = pd.Series(loops)
    loops2_sorted = loops2.sort_values()
    index_sorted = loops2_sorted.index.values
    resp = np.array(datas)[index_sorted]   
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

def create_loop_array_structured(data_case):
    resp = []
    for data in data_case:
        resp.append(data['loop_array'])
    
    return np.array(resp)

datas_case1 = organize_cases_by_loop(datas_case1, loops1)
datas_case2 = organize_cases_by_loop(datas_case2, loops2)

case1_structured_loop_array = create_loop_array_structured(datas_case1)
case2_structured_loop_array = create_loop_array_structured(datas_case2)

case1_oil_production = case1_structured_loop_array['oil_production'].flatten()
case2_oil_production = case2_structured_loop_array['oil_production'].flatten()

case1_gas_production = case1_structured_loop_array['gas_production'].flatten()
case2_gas_production = case2_structured_loop_array['gas_production'].flatten()

case1_oil_rate = case1_structured_loop_array['oil_rate'].flatten()
case2_oil_rate = case2_structured_loop_array['oil_rate'].flatten()

case1_gas_rate = case1_structured_loop_array['gas_rate'].flatten()
case2_gas_rate = case2_structured_loop_array['gas_rate'].flatten()

case1_simulation_time = case1_structured_loop_array['simulation_time'].flatten()
case2_simulation_time = case2_structured_loop_array['simulation_time'].flatten()

case1_time = case1_structured_loop_array['t'].flatten()
case2_time = case2_structured_loop_array['t'].flatten()
case1_time = case1_time/86400
case2_time = case2_time/86400


n_volumes_update = case2_structured_loop_array['n_volumes_update_base_functions'].flatten()
total_volumes_updated = case2_structured_loop_array['total_volumes_updated'].flatten()
case2_active_volumes = case2_structured_loop_array['active_volumes'].flatten()


max_time_case2 = case2_time.max()
max_time_case1 = case1_time.max()
max_time = min([case1_time.max(), case2_time.max()])

test1 = case1_time <= max_time
test2 = case2_time <= max_time

case1_time = case1_time[test1]
case1_gas_production = case1_gas_production[test1]
case1_oil_production = case1_oil_production[test1]
case1_simulation_time = case1_simulation_time[test1]
case1_oil_rate = case1_oil_rate[test1]
case1_gas_rate = case1_gas_rate[test1]

case2_time = case2_time[test2]
case2_gas_production = case2_gas_production[test2]
case2_oil_production = case2_oil_production[test2]
case2_simulation_time = case2_simulation_time[test2]
case2_oil_rate = case2_oil_rate[test2]
case2_gas_rate = case2_gas_rate[test2]

n_volumes_update = n_volumes_update[test2]
total_volumes_updated = total_volumes_updated[test2]
case2_active_volumes = case2_active_volumes[test2]

# import pdb; pdb.set_trace()


# import pdb; pdb.set_trace()


# p1 = extrair_dado(datas_case1, 'pressure')
# p2 = extrair_dado(datas_case2, 'pressure')

# erro = erro_abs(p1, p2)
# erro_rel = erro / p1

# erro_l2 = np.linalg.norm(erro_rel, axis=1)
# erro_max = np.max(erro_rel, axis=1)

########################################################
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
ax2.legend()

# plt.subplots_adjust(left=0.1,
#                     bottom=0.1, 
#                     right=0.9, 
#                     top=0.9, 
#                     wspace=0.4, 
#                     hspace=0.4)

fig.tight_layout()

plt.savefig(fig_str + description2 + 'Production' + '.png')
########################################################

######################################################
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
# plt.subplots_adjust(
#     top=0.9,
#     bottom=0.3,
#     wspace=0.4,
#     hspace=0.4
# )
    
# fig.tight_layout()
fig.suptitle('Volumes para atualizar as funcoes de base')

plt.savefig(fig_str + description2 + 'Total_volumes' + '.png')
####################################################

###################################################
fig, (ax1, ax2) = plt.subplots(2, 1)
ax1.plot(case2_time, case2_active_volumes, '-')
# ax1.fill(case2_time, case2_active_volumes)
# ax1.bar(case2_time, case2_active_volumes)
ax1.set_xlabel('time [days]')
ax1.set_ylabel('Active volumes')
ax1.set_xlim(min(case2_time), max(case2_time) + 1)

ax2.plot(case1_time, case1_simulation_time, '-', label='Finescale')
ax2.plot(case1_time, np.repeat(np.mean(case1_simulation_time), len(case1_simulation_time)), 0.05, color='black')
ax2.plot(case2_time, case2_simulation_time, '-', label='Adm')
ax2.plot(case2_time, np.repeat(np.mean(case2_simulation_time), len(case2_simulation_time)), 0.05, color='black')
# ax1.fill(case2_time, case2_active_volumes)
ax2.set_xlabel('time [days]')
ax2.set_ylabel('Simulation_time [s]')
ax2.set_xlim(min(case2_time), max(case2_time) + 1)
ax2.legend()
# fig.tight_layout()

plt.subplots_adjust(left=0.15,
                    bottom=0.1, 
                    right=0.9, 
                    top=0.9, 
                    wspace=0.4, 
                    hspace=0.4)

plt.savefig(fig_str + description2 + 'Active_volumes' + '.png')
#################################################

#################################################
fig, (ax1, ax2) = plt.subplots(2, 1)
ax1.plot(case1_time, case1_oil_rate, '-', label='Finescale')
ax1.plot(case2_time, case2_oil_rate, '-', label='Adm')
# ax1.fill(case2_time, case2_active_volumes)
# ax1.bar(case2_time, case2_active_volumes)
ax1.set_xlabel('time [days]')
ax1.set_ylabel('Oil rate [m3/s]')
ax1.set_xlim(min(case2_time), max(case2_time) + 1)
ax1.legend()

ax2.plot(case1_time, case1_gas_rate, '-', label='Finescale')
ax2.plot(case2_time, case2_gas_rate, '-', label='Adm')
ax2.set_xlim(min(case2_time), max(case2_time) + 1)
ax2.legend()
# fig.tight_layout()

plt.subplots_adjust(left=0.15,
                    bottom=0.1, 
                    right=0.9, 
                    top=0.9, 
                    wspace=0.4, 
                    hspace=0.4)

plt.savefig(fig_str + description2 + 'Flow_rate' + '.png')
###############################################



import pdb; pdb.set_trace()


        
    