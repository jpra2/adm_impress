import numpy as np
from scipy import io as sio
import os
import matplotlib.pyplot as plt

data_path = 'data'
# file_pressure_impress = os.path.join(data_path, 'pressures_p7.npy')
# file_pressure_impress = os.path.join(data_path, 'pressures_p7_2.npy')
file_pressure_impress = os.path.join(data_path, 'pressures_p7_3.npy')
# file_pressures_impress_g = os.path.join(data_path, 'pressures_p7_g.npy')
# file_pressures_impress_g = os.path.join(data_path, 'pressures_p7_g_2.npy')
file_pressures_impress_g = os.path.join(data_path, 'pressures_p7_g_3.npy')
file_centroids_imp = os.path.join(data_path, 'centroids.npy')
file_dados_mat = os.path.join(data_path, 'dados_p7.mat')
file_dados_mat_g = os.path.join(data_path, 'dados_p7_g.mat')
# file_sat_impress = os.path.join(data_path, 'saturations_p7.npy')
# file_sat_impress = os.path.join(data_path, 'saturations_p7_2.npy')
file_sat_impress = os.path.join(data_path, 'saturations_p7_3.npy')
# file_sat_g_impress = os.path.join(data_path, 'saturations_p7_g.npy')
# file_sat_g_impress = os.path.join(data_path, 'saturations_p7_g_2.npy')
file_sat_g_impress = os.path.join(data_path, 'saturations_p7_g_3.npy')

nl = 87

pressures_p7 = np.load(file_pressure_impress)[1:]
pressures_p7_g = np.load(file_pressures_impress_g)[1:]
sat_w_impress = np.load(file_sat_impress)
sat_w_g_impress = np.load(file_sat_g_impress)
pressures_p7 = pressures_p7[0:nl]
pressures_p7_g = pressures_p7_g[0:nl]
sat_w_impress = sat_w_impress[0:nl]
sat_w_g_impress = sat_w_g_impress[0:nl]
centroids_impress = np.load(file_centroids_imp)
dados_p7 = sio.loadmat(file_dados_mat)['dados']
dados_p7_g = sio.loadmat(file_dados_mat_g)['dados']


centroids_mat = dados_p7['centroids'][0][0]
centroids_mat[:,2] = 20 - centroids_mat[:,2]
pressures_p7_mat = dados_p7['pressure'][0]
pressures_p7_g_mat = dados_p7_g['pressure'][0]
sat_water_p7 = []
saturations_p7_mat = dados_p7['s'][0]

for sat in saturations_p7_mat:
    sat_water_p7.append(sat[:,0])

sat_water_p7 = np.array(sat_water_p7)

sat_water_p7_g = []
saturations_p7_g_mat = dados_p7_g['s'][0]

for sat in saturations_p7_g_mat:
    sat_water_p7_g.append(sat[:,0])

sat_water_p7_g = np.array(sat_water_p7_g)

def organize_data_by_centroids(centroids_ref, centroids_to_organize):

    n_centroids = len(centroids_ref)
    indexes = np.arange(n_centroids)
    new_index = []

    c2 = centroids_to_organize

    for i, j, k in zip(c2[:,0], c2[:,1], c2[:,2]):
        v1 = centroids_ref[:,0] == i
        v2 = centroids_ref[:,1] == j
        v3 = centroids_ref[:,2] == k
        vf = v1 & v2 & v3
        new_index.append(indexes[vf])

    new_index = np.array(new_index).flatten()
    return new_index

organize_index = organize_data_by_centroids(centroids_impress, centroids_mat)

def organize_data_list_by_index(data, index):

    new_data = []

    for line in data:
        new_data.append(line[index].flatten()[::-1])

    return np.array(new_data)

pressures_p7_mat_comp = organize_data_list_by_index(pressures_p7_mat, organize_index)
pressures_p7_g_mat_comp = organize_data_list_by_index(pressures_p7_g_mat, organize_index)
sat_water_p7_mat_comp = organize_data_list_by_index(sat_water_p7, organize_index)
sat_water_p7_g_mat_comp = organize_data_list_by_index(sat_water_p7_g, organize_index)

pressures_p7_mat_comp = pressures_p7_mat_comp[0:nl]
pressures_p7_g_mat_comp = pressures_p7_g_mat_comp[0:nl]
sat_water_p7_mat_comp = sat_water_p7_mat_comp[0:nl]
sat_water_p7_g_mat_comp = sat_water_p7_g_mat_comp[0:nl]

def get_erros(data_comp, correct_data):

    linf = []
    l2 = []

    for d1, d2 in zip(data_comp, correct_data):
        erro_abs = get_abs_err(d1, d2)
        ninf = get_ninf(erro_abs, d2)
        nl2 = get_nl2(erro_abs, d2)
        linf.append(ninf)
        l2.append(nl2)
    
    return np.array(linf), np.array(l2)

def get_abs_err(dat1, dat2):
    return np.absolute(dat1 - dat2)

def get_ninf(err_abs, correct_data):
    return (err_abs / (max(np.absolute(correct_data)))).max()

def get_nl2(err_abs, correct_data):
    return (np.linalg.norm(err_abs))/(np.linalg.norm(correct_data))

linf_p7, l2_p7 = get_erros(pressures_p7, pressures_p7_mat_comp)
linf_p7_g, l2_p7_g = get_erros(pressures_p7_g, pressures_p7_g_mat_comp)

linf_sat_p7, l2_sat_p7 = get_erros(sat_w_impress, sat_water_p7_mat_comp)
linf_sat_p7[0] = 0
l2_sat_p7[0] = 0
linf_sat_p7_g, l2_sat_p7_g = get_erros(sat_w_g_impress, sat_water_p7_g_mat_comp)
linf_sat_p7_g[0] = 0
l2_sat_p7_g[0] = 0

nums = np.arange(len(linf_p7))

fig, ax = plt.subplots(nrows=2, ncols=2)
# fig.tight_layout()
plt.subplots_adjust(hspace=0.5, wspace=0.4)

# ax[0][0] = fig.add_subplot(2,2,1)
ax[0][0].plot(nums, linf_p7)
ax[0][0].set_title('Pressão')
ax[0][0].set_xlabel('Passos de tempo')
ax[0][0].set_ylabel('$L_{\infty}$', fontsize=12)
ax[0][0].set_ylim([0, linf_p7.max()*(1.1)])
ax[0][0].set_xlim([0, nums.max()*(1.1)])

# ax[0][1] = fig.add_subplot(2,2,2)
ax[0][1].plot(nums, l2_p7)
ax[0][1].set_title('Pressão')
ax[0][1].set_xlabel('Passos de tempo')
ax[0][1].set_ylabel('$L_{2}$', fontsize=12)
ax[0][1].set_ylim([0, l2_p7.max()*(1.1)])
ax[0][1].set_xlim([0, nums.max()*(1.1)])

# ax[1][0] = fig.add_subplot(2,2,3)
ax[1][0].plot(nums, linf_sat_p7)
ax[1][0].set_title('Saturação')
ax[1][0].set_xlabel('Passos de tempo')
ax[1][0].set_ylabel('$L_{\infty}$', fontsize=12)
ax[1][0].set_ylim([0, linf_sat_p7.max()*(1.1)])
ax[1][0].set_xlim([0, nums.max()*(1.1)])

# ax[1][1] = fig.add_subplot(2,2,4)
ax[1][1].plot(nums, l2_sat_p7)
ax[1][1].set_title('Saturaçao')
ax[1][1].set_xlabel('Passos de tempo')
ax[1][1].set_ylabel('$L_{2}$', fontsize=12)
ax[1][1].set_ylim([0, l2_sat_p7.max()*(1.1)])
ax[1][1].set_xlim([0, nums.max()*(1.1)])

plt.suptitle('Vazão pressão sem gravidade')
plt.savefig(os.path.join('data', 'graf3.png'))
# plt.show()

fig, ax = plt.subplots(nrows=2, ncols=2)
# fig.tight_layout()
plt.subplots_adjust(hspace=0.5, wspace=0.4)

ax[0][0].plot(nums, linf_p7_g)
ax[0][0].set_title('Pressão')
ax[0][0].set_xlabel('Passos de tempo')
ax[0][0].set_ylabel('$L_{\infty}$', fontsize=12)
ax[0][0].set_ylim([0, linf_p7_g.max()*(1.1)])
ax[0][0].set_xlim([0, nums.max()*(1.1)])

ax[0][1].plot(nums, l2_p7_g)
ax[0][1].set_title('Pressão')
ax[0][1].set_xlabel('Passos de tempo')
ax[0][1].set_ylabel('$L_{2}$', fontsize=12)
ax[0][1].set_ylim([0, l2_p7_g.max()*(1.1)])
ax[0][1].set_xlim([0, nums.max()*(1.1)])

ax[1][0].plot(nums, linf_sat_p7_g)
ax[1][0].set_title('Saturação')
ax[1][0].set_xlabel('Passos de tempo')
ax[1][0].set_ylabel('$L_{\infty}$', fontsize=12)
ax[1][0].set_ylim([0, linf_sat_p7_g.max()*(1.1)])
ax[1][0].set_xlim([0, nums.max()*(1.1)])

ax[1][1].plot(nums, l2_sat_p7_g)
ax[1][1].set_title('Saturação')
ax[1][1].set_xlabel('Passos de tempo')
ax[1][1].set_ylabel('$L_{2}$', fontsize=12)
ax[1][1].set_ylim([0, l2_sat_p7_g.max()*(1.1)])
ax[1][1].set_xlim([0, nums.max()*(1.1)])

plt.suptitle('Vazão pressão com gravidade')
plt.savefig(os.path.join('data', 'graf4.png'))
# plt.show()









# import pdb; pdb.set_trace()