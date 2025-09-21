from packs.manager import SimulationData
import matplotlib.pyplot as plt
from matplotlib import markers
import matplotlib.colors as mcolors
from typing import Sequence
from packs import defpaths
import os
import numpy as np

from packs.multiscale.unstructured.test.test_brazil import (  
    f_l2_error_paper_artur,
    f_linf_percent_paper_artur, 
    f_l2_error_percent_paper_artur
)

from packs.examples.same_functions import (
    load_ps_results
)


def get_markers() -> list:
    # markers = {'': 'nothing', '*': 'star', '+': 'plus', ',': 'pixel', '.': 'point', '1': 'tri_down', '2': 'tri_up', '3': 'tri_left', '4': 'tri_right', '8': 'octagon', '<': 'triangle_left', '>': 'triangle_right', 'D': 'diamond', 'H': 'hexagon2', 'P': 'plus_filled', 'X': 'x_filled', '^': 'triangle_up', '_': 'hline', 'd': 'thin_diamond', 'h': 'hexagon1', 'none': 'nothing', 'o': 'circle', 'p': 'pentagon', 's': 'square', 'v': 'triangle_down', 'x': 'x', '|': 'vline', 0: 'tickleft', 1: 'tickright', 10: 'caretupbase', 11: 'caretdownbase', 2: 'tickup', 3: 'tickdown', 4: 'caretleft', 5: 'caretright', 6: 'caretup', 7: 'caretdown', 8: 'caretleftbase', 9: 'caretrightbase'}
    markers = {'+': 'plus', '.': 'point', '1': 'tri_down', '2': 'tri_up', '3': 'tri_left', '4': 'tri_right', '8': 'octagon', '<': 'triangle_left', '>': 'triangle_right', 'D': 'diamond', 'H': 'hexagon2', 'P': 'plus_filled', 'X': 'x_filled', '^': 'triangle_up', '_': 'hline', 'd': 'thin_diamond', 'h': 'hexagon1', 'none': 'nothing', 'o': 'circle', 'p': 'pentagon', 's': 'square', 'v': 'triangle_down', 'x': 'x', '|': 'vline', 0: 'tickleft', 1: 'tickright', 10: 'caretupbase', 11: 'caretdownbase', 2: 'tickup', 3: 'tickdown', 4: 'caretleft', 5: 'caretright', 6: 'caretup', 7: 'caretdown', 8: 'caretleftbase', 9: 'caretrightbase'}
    return list(markers.keys())

def get_linestyle_str() -> list:
    linestyle_str = [
     ('solid', 'solid'),      # Same as (0, ()) or '-'
     ('dotted', 'dotted'),    # Same as ':'
     ('dashed', 'dashed'),    # Same as '--'
     ('dashdot', 'dashdot')
    ]  # Same as '-.'

    linestyles = []
    for line in linestyle_str:
        linestyles.append(line[0])
    
    return linestyles

def get_linestyle_tuple() -> list:
    linestyle_tuple = [
     ('loosely dotted',        (0, (1, 10))),
     ('dotted',                (0, (1, 5))),
     ('densely dotted',        (0, (1, 1))),

     ('long dash with offset', (5, (10, 3))),
     ('loosely dashed',        (0, (5, 10))),
     ('dashed',                (0, (5, 5))),
     ('densely dashed',        (0, (5, 1))),

     ('loosely dashdotted',    (0, (3, 10, 1, 10))),
    #  ('dashdotted',            (0, (3, 5, 1, 5))),
     ('dashdotted',            (0, (5, 3, 1, 3))),
     ('densely dashdotted',    (0, (3, 1, 1, 1))),

     ('dashdotdotted',         (0, (3, 5, 1, 5, 1, 5))),
     ('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10))),
     ('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1)))]
    
    linestyles = []
    for line in linestyle_tuple:
        linestyles.append(line[1])
    
    return linestyles

def basic_colors() -> list:
    return list(mcolors.BASE_COLORS.keys())

def tableu_colors() -> list:
    return list(mcolors.TABLEAU_COLORS.keys())

def css4_colors() -> list:
    return list(mcolors.CSS4_COLORS.keys())

def xkcd_colors() -> list:
    return list(mcolors.XKCD_COLORS.keys())


def plot_cum_oil(finescale_sim: SimulationData, nuadm_sims: Sequence[SimulationData], figname: str):

    x_label = 'PVI'
    y_label = 'Cumulative Oil'
    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot()

    """
    malha fina : continuo
    nuadm: dotted
    marcador distingue as simulacoes nuadm
    """
    markers = get_markers()
    linestyle_nuadm = get_linestyle_str()[1]
    basic_colors_str = basic_colors()
    basic_colors_str.remove('k')
    n_nuadmsims = len(nuadm_sims)
    markers_size = np.arange(7, 7+n_nuadmsims)

    nplus = 2
    for i in range(1, markers_size.shape[0]):
        markers_size[i] += nplus
    
    markers_size = markers_size[::-1]

    
    for i, sim in enumerate(nuadm_sims):
        max_loop = sim['all_loops'].max()
        cum_oil = -1*sim['all_cumulative_oil']
        vpi = sim['all_vpi']
        color = basic_colors_str[i]
        ax.plot(
            vpi, 
            cum_oil, 
            # label=sim.label, 
            label='NU-ADM', 
            # marker=markers[i],
            # markersize=markers_size[i],
            linewidth=4, 
            # linestyle=linestyle_nuadm,
            color=color)
    
    test_finescale = finescale_sim['all_loops'] <= max_loop
    
    ax.plot(
        finescale_sim['all_vpi'][test_finescale],
        -1*finescale_sim['all_cumulative_oil'][test_finescale],
        label='finescale',
        color='k'
    )

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.legend(handlelength=5)

    path_fig = os.path.join(defpaths.plots_folder, figname)
    fig.savefig(path_fig, dpi=500)
    
def plot_wor(finescale_sim: SimulationData,  nuadm_sims: Sequence[SimulationData], figname: str):
    x_label = 'PVI'
    y_label = 'WOR'
    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot()

    """
    malha fina : continuo
    nuadm: dotted
    marcador distingue as simulacoes nuadm
    """
    markers = get_markers()
    linestyle_nuadm = get_linestyle_str()[1]
    basic_colors_str = basic_colors()
    basic_colors_str.remove('k')
    n_nuadmsims = len(nuadm_sims)
    markers_size = np.arange(7, 7+n_nuadmsims)

    nplus = 2
    for i in range(1, markers_size.shape[0]):
        markers_size[i] += nplus
    
    markers_size = markers_size[::-1]

    
    for i, sim in enumerate(nuadm_sims):
        max_loop = sim['all_loops'].max()
        water_flux = sim['water_flux']
        oil_flux = sim['oil_flux']
        wor = water_flux/oil_flux
        vpi = sim['all_vpi']
        color = basic_colors_str[i]
        ax.plot(
            vpi, 
            wor, 
            # label=sim.label, 
            label='NU-ADM', 
            # marker=markers[i+1],
            # markersize=markers_size[i], 
            # linestyle=linestyle_nuadm,
            linewidth=3,
            color=color)
    
    test_finescale = finescale_sim['all_loops'] <= max_loop
    wor_finescale = finescale_sim['water_flux']/finescale_sim['oil_flux']
    
    ax.plot(
        finescale_sim['all_vpi'][test_finescale],
        wor_finescale[test_finescale],
        label='finescale',
        color='k'
    )

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.legend(handlelength=5)

    path_fig = os.path.join(defpaths.plots_folder, figname)
    fig.savefig(path_fig, dpi=500)

def plot_nuadm_percent(finescale_sim: SimulationData,  nuadm_sims: Sequence[SimulationData], figname: str):
    x_label = 'PVI'
    y_label = '%NUADM VOLUMES'
    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot()

    """
    malha fina : continuo
    nuadm: dotted
    marcador distingue as simulacoes nuadm
    """
    markers = get_markers()
    linestyle_nuadm = get_linestyle_str()[1]
    basic_colors_str = basic_colors()
    basic_colors_str.remove('k')
    n_nuadmsims = len(nuadm_sims)
    markers_size = np.arange(7, 7+n_nuadmsims)

    nplus = 2
    for i in range(1, markers_size.shape[0]):
        markers_size[i] += nplus
    
    markers_size = markers_size[::-1]

    pressure0_finescale = finescale_sim['pressure_0']
    nfinevolumes = pressure0_finescale.shape[0]

    
    for i, sim in enumerate(nuadm_sims):
        nuadm_volumes = 100*sim['nuadm_vols']/nfinevolumes
        vpi = sim['all_vpi']
        color = basic_colors_str[i]
        ax.plot(
            vpi, 
            nuadm_volumes, 
            # label=sim.label, 
            # marker=markers[i+1],
            # markersize=markers_size[i], 
            # linestyle=linestyle_nuadm,
            linewidth=3,
            color=color)

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    # ax.legend(handlelength=5)

    path_fig = os.path.join(defpaths.plots_folder, figname)
    fig.savefig(path_fig, dpi=500)


def chueh_test():

    x = np.linspace(0, 1, 800)
    y = 0.5 + 0.1*np.sin(10*x)
    y2 = y + 0.12
    y3 = y - 0.12

    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.plot(x, y)
    ax.plot(x, y2)
    ax.plot(x, y3)
    fig_path = os.path.join(defpaths.plots_folder, 'chueh_1.svg')
    fig.savefig(fig_path)


def plot_errors(vpis_to_plot, finescale_sim, nuadm_sim):
    max_delta = 1e-9
    fig_sat_name = 'saturation_error.png'
    fig_pressure_name = 'pressure_error.png'
    
    pres_str = 'pressure_'
    sat_str = 'saturation_'
    
    pressure_linf = []
    pressure_l2 = []
    
    sat_linf = []
    sat_l2 = []
    
    all_vpis_finescale = finescale_sim['all_vpi']
    all_vpi_nuadm = nuadm_sim['all_vpi']
    
    all_loops_finescale = finescale_sim['all_loops']
    all_loops_coarse = nuadm_sim['all_loops']
    
    fs = finescale_sim
    cs = nuadm_sim


    # vpi_test = 0.04
    # vpi_test = 0.12
    # vpi_test = 0.24
    # test1 = np.absolute(all_vpi_nuadm - vpi_test) <= max_delta
    # test2 = np.absolute(all_vpis_finescale - vpi_test) <= max_delta
    # loop_test = all_loops_coarse[test1]
    # loop_finescale = all_loops_finescale[test2]
    
    # # loop_t = 981
    # # vpi = all_vpi_nuadm[all_loops_coarse==loop_t]
    
    # max_loop = min([all_loops_finescale.max(), all_loops_coarse.max()])
    
    marki = -1
    
    for i, vpi in enumerate(vpis_to_plot):
        test_fs = np.absolute(all_vpis_finescale - vpi) <= max_delta
        test_nuadm = np.absolute(all_vpi_nuadm - vpi) <= max_delta
        marki = i
        
        
        if np.any(test_fs) and np.any(test_nuadm):
            pass
        else:
            break
            
        
        loop_fs = finescale_sim['all_loops'][test_fs][0]
        loop_nuadm = nuadm_sim['all_loops'][test_nuadm][0]
        
        # key_pressure_fs = pres_str + str(loop_fs)
        # key_pressure_nuadm = pres_str + str(loop_nuadm)
        
        # key_sat_fs = sat_str + str(loop_fs)
        # key_sat_nuadm = sat_str + str(loop_nuadm)
        
        # pressure_fs = finescale_sim[key_pressure_fs]
        # pressure_nuadm = nuadm_sim[key_pressure_nuadm]
        
        pressure_fs, saturation_fs = load_ps_results(loop_fs, defpaths.pressure_results, defpaths.saturation_results)
        pressure_nuadm, saturation_nuadm = load_ps_results(loop_nuadm, defpaths.pressure_results_ms, defpaths.saturation_results_ms)
        
        # saturation_fs = finescale_sim[key_sat_fs]
        # saturation_nuadm = nuadm_sim[key_sat_nuadm]
        
        linf_pressure = f_linf_percent_paper_artur(pressure_fs, pressure_nuadm)/100
        err2_pressure = f_l2_error_paper_artur(pressure_fs, pressure_nuadm)
        
        linf_sat = f_linf_percent_paper_artur(saturation_fs, saturation_nuadm)/100
        err2_sat = f_l2_error_paper_artur(saturation_fs, saturation_nuadm)
        
        pressure_linf.append(linf_pressure)
        pressure_l2.append(err2_pressure)
        
        sat_linf.append(linf_sat)
        sat_l2.append(err2_sat)
    
    pressure_linf = np.array(pressure_linf)
    pressure_l2 = np.array(pressure_l2)
    
    sat_linf = np.array(sat_linf)[1:]
    sat_l2 = np.array(sat_l2)[1:]
    
    markers = get_markers()
    linestyles = get_linestyle_str()
    basic_colors_str = basic_colors()
    # basic_colors_str.remove('k')
    # n_nuadmsims = len(nuadm_sims)
    n_nuadmsims = 1
    markers_size = np.arange(7, 7+n_nuadmsims)

    nplus = 2
    for i in range(1, markers_size.shape[0]):
        markers_size[i] += nplus
    
    markers_size = markers_size[::-1]
    plt.clf()
    plt.rcParams['text.usetex'] = True
    fig = plt.figure()
    ax = fig.add_subplot()
    
    ax.plot(
        vpis_to_plot[0:marki],
        pressure_linf,
        label=r'$||p||_{\infty}$', 
        marker=markers[0],
        markersize=markers_size[0], 
        color='k'
    )
    ax.plot(
        vpis_to_plot[0:marki],
        pressure_l2,
        label=r'$||p||_{2}$', 
        marker=markers[1],
        markersize=markers_size[0], 
        color='red'
    )    
    ax.set_yscale('log')
    
    desired_yticks = [0.01, 0.1, 1]
    ax.set_yticks(desired_yticks)
    
    ax.set_xlabel('PVI')
    # ax.set_ylabel(r'')
    ax.legend(handlelength=5)

    path_fig = os.path.join(defpaths.plots_folder, 'pressure_error_biphasic.png')
    fig.savefig(path_fig, dpi=500)
    
    plt.clf()
    # plt.rcParams['text.usetex'] = True
    fig = plt.figure()
    ax = fig.add_subplot()
    
    ax.plot(
        vpis_to_plot[1:marki],
        sat_linf,
        label=r'$||S_{w}||_{\infty}$', 
        marker=markers[0],
        markersize=markers_size[0], 
        color='k'
    )
    ax.plot(
        vpis_to_plot[1:marki],
        sat_l2,
        label=r'$||S_{w}||_{2}$', 
        marker=markers[1],
        markersize=markers_size[0], 
        color='red'
    )    
    ax.set_yscale('log')
    
    desired_yticks = [0.01, 0.1, 1]
    ax.set_yticks(desired_yticks)
    
    ax.set_xlabel('PVI')
    # ax.set_ylabel(r'')
    ax.legend(handlelength=5)

    path_fig = os.path.join(defpaths.plots_folder, 'saturation_error_biphasic.png')
    fig.savefig(path_fig, dpi=500)
    
    
    
    
        
    import pdb; pdb.set_trace()
        
        
def plot_err_layers():
    err = np.load('err.npy')
    it = np.load('it.npy')
    
    err2 = np.load('err2.npy')
    it2 = np.load('it2.npy')
    
    markers = get_markers()
    linestyles = get_linestyle_str()
    basic_colors_str = basic_colors()
    # basic_colors_str.remove('k')
    # n_nuadmsims = len(nuadm_sims)
    
    plt.clf()
    plt.rcParams['text.usetex'] = True
    plt.rcParams['legend.fontsize'] = 16
    plt.rcParams['axes.labelsize'] = 16
    plt.rcParams["font.family"] = "Times New Roman"
    fig = plt.figure()
    ax = fig.add_subplot()
    
    ax.plot(
        it,
        err,
        label='AMS-U', 
        marker=markers[1], 
        color='k'
    )
    ax.plot(
        it2,
        err2,
        label='NU-ADM-U', 
        marker=markers[1], 
        color='red'
    )    
    ax.set_yscale('log')
    ax.set_ylim(1e-13, 10)
    
    ax.set_xlabel('Number of iterations')
    ax.set_ylabel(r'$||\bf{r}||_{2}$')
    ax.legend(handlelength=5)
    path_fig = os.path.join(defpaths.plots_folder, 'layers_conv.png')
    fig.savefig(path_fig, dpi=500)
    
    
        
        
        
        
        
    
    
    










    
    


def plot_graphs():

    finescale_sim = SimulationData('biphasic_ameba_finescale4')
    # finescale_sim = SimulationData('biphasic_het1_finescale')
    # finescale_sim = SimulationData('biphasic_het_coarse1')
    # finescale_sim = SimulationData('biphasic_sin_chueh_fine1')
    finescale_sim.load_data()
    fs = finescale_sim
    # vpis_to_plot = np.linspace(0, 0.24, 25)
    # vpis_to_plot = np.linspace(0, 0.1155, 16)
    vpis_to_plot = np.linspace(0, 0.6, 31)

    nuadm_sims_str = ['biphasic_ameba_coarse4']
    # nuadm_sims_str = ['biphasic_het_coarse1_1']
    # nuadm_sims_str = ['biphasic_sin1_coarse2']
    nuadm_sims = []
    for i, name in enumerate(nuadm_sims_str):
        data_sim = SimulationData(name)
        data_sim.load_data()
        if data_sim.verify_name_in_data_names('label') == True:
            pass
        else:
            data_sim.insert_or_update_data({'label': np.array(['nuadm'+str(i)])})
        nuadm_sims.append(data_sim)
    
    
    
    # plot_cum_oil(finescale_sim, nuadm_sims, 'cum_oil_ameba_new.png')
    # plot_wor(finescale_sim, nuadm_sims, 'wor_ameba_new.png')
    # plot_nuadm_percent(finescale_sim, nuadm_sims, 'nuadm_percent_ameba_new.png')
    # plot_errors(vpis_to_plot, finescale_sim, nuadm_sims[0])
    
    plot_err_layers()
    
    import pdb; pdb.set_trace()





    pass
