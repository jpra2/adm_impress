from packs.data_class.compositional_cumulative_datamanager import CumulativeCompositionalDataManager
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from packs.cases.compositional_adm_cases.compressible_oil.all_functions import organize_cases_by_loop, extrair_dado, erro_abs, get_data_from_loop_array, create_loop_array_structured
from packs.cases.compositional_adm_cases.compressible_oil import descriptions
import os

# key_list = ['loop_array', 'pressure']
key_list = ['loop_array']

# description = 'case1_finescale_'
# description1 = 'case3_finescale_3k'
# description1 = 'case17_finescale_6k_5000_'
# description1 = 'case23_finescale_6k_5000_' # iterative
# description1 = 'case37_finescale_80x80_Firoozabadi_direct_solver' # iterative
description1 = descriptions.case1_finescale # iterative
# description = 'case2_adm_'
# compositional_data = CompositionalData(description=description)
case1 = CumulativeCompositionalDataManager(description=description1)
datas_case1 = case1.load_all_datas_from_keys(key_list)

fig_str = 'figura_'
# fig_ext = '.svg'
fig_ext = '.png'

# description = 'case2_adm_'
# description = 'case4_adm_3k'
# description = 'case5_adm_3k'
# description = 'case7_adm_3k'
# description = 'case8_adm_3k'
# description2 = 'case9_adm_3k'
# description2 = 'case10_adm_3k'
# description2 = 'case11_adm_3k'
# description2 = 'case18_adm_6k_5000_'
# description2 = 'case19_adm_6k_5000_'
# description2 = 'case20_adm_6k_5000_'
# description2 = 'case21_adm_6k_5000_'
# description2 = 'case36_adm_80x80_Firoo_tams_solver_new_prolong'
# description2 = 'case37_adm_80x80_Firoo_tams_solver_new_prolong_coarsewells_level0'
# description2 = 'case38_adm_80x80_Firoo_iterative_CG_new_prolong_coarsewells_level0'
# description2 = 'case39_adm_80x80_Firoo_iterative_CG_new_prolong_coarsewells_level0_cr-10'
# description2 = 'case45_test'
# description2 = 'case46_test-limit_maxiter_tams-10'
# description2 = 'case47-test-limit_maxiter_tams_cg-20_thres-update-BF-0.3'
# description2 = 'case48-test-limit_maxiter_tams_cg-1000_thres-update-BF-0.1_by_dvtol'
# description2 = 'case49-test-limit_maxiter_tams_cg-1000_thres-update-BF-0.1_by_dvtol'
description2 = descriptions.case2_adm_description
# description2 = 'case22_adm_6k_5000_'
# description2 = 'case23_finescale_6k_5000_'
case2 = CumulativeCompositionalDataManager(description=description2)
# datas_case2 = case2.load_all_datas()
datas_case2 = case2.load_all_datas_from_keys(key_list)

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
# max_time = 0.7

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
plt.clf()
fig, ax1 = plt.subplots()

ax1.plot(case1_time, case1_oil_production, '-', label='Finescale')
ax1.plot(case2_time, case2_oil_production, '-', label='Adm')
ax1.set_ylabel('Oil production[m³]')
ax1.set_xlabel('time [days]')
ax1.legend()
plt.savefig(fig_str + description2 + 'Oil_production' + fig_ext)

plt.clf()
fig, ax1 = plt.subplots()

ax1.plot(case1_time, case1_gas_production, '-', label='Finescale')
ax1.plot(case2_time, case2_gas_production, '-', label='Adm')
ax1.set_ylabel('Gas production[m³]')
ax1.set_xlabel('time [days]')
ax1.legend()
plt.savefig(fig_str + description2 + 'Gas_production' + fig_ext)
# plt.subplots_adjust(left=0.1,
#                     bottom=0.1,
#                     right=0.9,
#                     top=0.9,
#                     wspace=0.4,
#                     hspace=0.4)

fig.tight_layout()

# plt.savefig(fig_str + description2 + 'Production' + '.png')
#plt.savefig(fig_str + description2 + 'Production_iterative' + fig_ext)
########################################################

######################################################
plt.clf()
fig, ax1 = plt.subplots()
ax1.plot(case2_time, n_volumes_update, '-')
ax1.set_xlabel('time [days]')
ax1.set_ylabel('N volumes')

plt.savefig(fig_str + description2 + 'N_volumes' + fig_ext)

plt.clf()
fig, ax1 = plt.subplots()
ax1.plot(case2_time, total_volumes_updated, '-')
ax1.set_xlabel('time [days]')
ax1.set_ylabel('Total volumes')

plt.savefig(fig_str + description2 + 'Total volumes' + fig_ext)

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
#fig.suptitle('Volumes para atualizar as funcoes de base')

# plt.savefig(fig_str + description2 + 'Total_volumes' + '.png')
####################################################

###################################################
plt.clf()
fig, ax1 = plt.subplots()
ax1.plot(case2_time, case2_active_volumes, '-')
# ax1.fill(case2_time, case2_active_volumes)
# ax1.bar(case2_time, case2_active_volumes)
ax1.set_xlabel('time [days]')
ax1.set_ylabel('Active volumes')
ax1.set_xlim(min(case2_time), max(case2_time))
plt.savefig(fig_str + description2 + 'Active_volumes' + fig_ext)

plt.clf()
fig, ax1 = plt.subplots()
ax1.plot(case1_time, case1_simulation_time, '-', label='Finescale')
ax1.plot(case2_time, case2_simulation_time, '-', label='NU-ADM')
# ax2.plot(case1_time, np.repeat(np.mean(case1_simulation_time), len(case1_simulation_time)), 0.05, color='red')
# ax2.plot(case2_time, np.repeat(np.mean(case2_simulation_time), len(case2_simulation_time)), 0.05, color='black')
ax1.axhline(y=np.mean(case1_simulation_time), color='red', label='Finescale mean')
ax1.axhline(y=np.mean(case2_simulation_time), color='black', label='NU-ADM mean')
ax1.annotate(round(np.mean(case1_simulation_time), 1), (0, np.mean(case1_simulation_time)))
ax1.annotate(round(np.mean(case2_simulation_time), 1), (0, np.mean(case2_simulation_time)))
ax1.grid(True)
t1 = np.mean(np.mean(case1_simulation_time))
t2 = np.mean(case2_simulation_time)
t1sum = case1_simulation_time.sum()
t2sum = case2_simulation_time.sum()
# import pdb; pdb.set_trace()
# ax1.fill(case2_time, case2_active_volumes)
ax1.set_xlabel('time [days]')
ax1.set_ylabel('Simulation_time [s]')
ax1.set_xlim(min(case2_time), max(case2_time))
ax1.legend()
# fig.tight_layout()

plt.subplots_adjust(left=0.15,
                    bottom=0.1,
                    right=0.9,
                    top=0.9,
                    wspace=0.4,
                    hspace=0.4)
plt.savefig(fig_str + description2 + 'Tempo desimulação' 'Active_volumes' + '.png')

#################################################

#################################################
plt.clf()
fig, ax1 = plt.subplots()
ax1.plot(case1_time, case1_oil_rate, '-', label='Finescale')
ax1.plot(case2_time, case2_oil_rate, '-', label='Adm')
# ax1.fill(case2_time, case2_active_volumes)
# ax1.bar(case2_time, case2_active_volumes)
ax1.set_xlabel('time [days]')
ax1.set_ylabel('Oil rate [m3/s]')
ax1.set_xlim(min(case2_time), max(case2_time))
ax1.legend()
plt.savefig(fig_str + description2 + 'Oil_rate' + fig_ext)

fig, ax1 = plt.subplots()
ax1.plot(case1_time, case1_gas_rate, '-', label='Finescale')
ax1.plot(case2_time, case2_gas_rate, '-', label='Adm')
ax1.set_xlabel('time [days]')
ax1.set_ylabel('Gas rate [m3/s]')
ax1.set_xlim(min(case2_time), max(case2_time))
ax1.legend()
# fig.tight_layout()
plt.savefig(fig_str + description2 + 'Gas_rate' + fig_ext)
plt.subplots_adjust(left=0.15,
                    bottom=0.1,
                    right=0.9,
                    top=0.9,
                    wspace=0.4,
                    hspace=0.4)

# plt.savefig(fig_str + description2 + 'Flow_rate' + '.png')

###############################################

########################
case2_tams_iterations = case2_structured_loop_array['tams_iterations'].flatten()
case2_tams_iterations = case2_tams_iterations[test2]
plt.clf()
fig, ax = plt.subplots()
ax.plot(case2_time, case2_tams_iterations, '-', label='Tams Iterations')
ax.set_xlabel('time [days]')
ax.set_ylabel('Iterations')
# fig.suptitle('Volumes para atualizar as funcoes de base')
plt.savefig(fig_str + description2 + 'Tams_iterations' + fig_ext)
########################

#########################################
# fig, ax = plt.subplots()
plt.clf()
y_values = [int(case1_simulation_time.sum()), int(case2_simulation_time.sum())]
x_values = ['Finescale', 'NU-ADM']
bars = plt.bar(x_values, y_values, width=0.3)
plt.ylabel('Total simulation time')
for bar in bars:
    yval = bar.get_height()
    plt.text(bar.get_x(), yval + 100, yval)
plt.savefig(fig_str + description2 + 'Total_simulation_time' + fig_ext)
###############################

# assign your bars to a variable so their attributes can be accessed
# bars = plt.bar(x, height=y, width=.4)

# # access the bar attributes to place the text in the appropriate location
# for bar in bars:
#     yval = bar.get_height()
#     plt.text(bar.get_x(), yval + .005, yval)

###########################################
# # norma do erro production interpolado
# x_finescale = case1_time
# y_finescale = case1_oil_production
# x_adm = case2_time
# y_adm = case2_oil_production

def get_minmax_interp(v0, v1):
    n_points = min([len(v0), len(v1)])
    minx_interp = max([v0.min(), v1.min()])
    maxx_interp = min(v0.max(), v1.max())
    resp_interp = np.linspace(minx_interp, maxx_interp, n_points)
    return resp_interp 

def erros_interp(x_finescale, y_finescale, x_adm, y_adm, name=''):
    
    x_interp = get_minmax_interp(x_finescale, x_adm)
    
    y_interp_finescale = np.interp(x_interp, x_finescale, y_finescale)
    y_interp_adm = np.interp(x_interp, x_adm, y_adm)
    
    erro_abs = np.absolute(y_interp_finescale - y_interp_adm)
    erro_rel = erro_abs/y_interp_finescale
    log10_erro_rel = np.log10(erro_rel)
    
    plt.clf()
    # fig, ax1 = plt.subplots()
    # ax1.plot(x_interp, log10_erro_rel, '-', label='Log10 do erro relativo')
    plt.plot(x_interp, log10_erro_rel, '-')
    plt.ylabel('Log10 do erro relativo')
    plt.xlabel('time [days]')
    # ax1.set_ylabel('Log10 do erro relativo')
    # ax1.set_xlabel('time [days]')
    # ax1.legend()
    
    plt.savefig(fig_str + description2 + 'log10_erro_rel_' + name + fig_ext)
    plt.clf()
    return {
        'erro_abs': erro_abs,
        'erro_rel': erro_rel,
        'log10_erro_rel': log10_erro_rel,
        'x': x_interp
    }

file_name = 'erros_oil_prod_' + description2 + '.npz'
file_name = os.path.join('flying', file_name)
erros_oil = erros_interp(case1_time, case1_oil_production, case2_time, case2_oil_production, name='oil')
np.savez(file_name, **erros_oil)


file_name = 'erros_gas_prod_' + description2 + '.npz'
file_name = os.path.join('flying', file_name)
erros_gas = erros_interp(case1_time, case1_gas_production, case2_time, case2_gas_production, name='gas')
np.savez(file_name, **erros_gas)
########################




import pdb; pdb.set_trace()
