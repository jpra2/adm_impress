from packs.manager import SimulationData
import matplotlib.pyplot as plt
from matplotlib import markers
import matplotlib.colors as mcolors
from typing import Sequence
from packs import defpaths
import os
import numpy as np

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

    x_label = 'vpi'
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
        cum_oil = -1*sim['all_cumulative_oil']
        vpi = sim['all_vpi']
        color = basic_colors_str[i]
        ax.plot(
            vpi, 
            cum_oil, 
            label=sim.label, 
            marker=markers[i],
            markersize=markers_size[i], 
            linestyle=linestyle_nuadm,
            color=color)
    
    ax.plot(
        finescale_sim['all_vpi'],
        -1*finescale_sim['all_cumulative_oil'],
        label='finescale',
        color='k'
    )

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.legend(handlelength=5)

    path_fig = os.path.join(defpaths.plots_folder, figname)
    fig.savefig(path_fig)
    
    
    


    









    
    


def plot_graphs():

    finescale_sim = SimulationData('biphasic_layer_finescale4')
    finescale_sim.load_data()

    nuadm_sims_str = ['biphasic_layers_coarse4', 'biphasic_ameba_coarse4f']
    nuadm_sims = []
    for i, name in enumerate(nuadm_sims_str):
        data_sim = SimulationData(name)
        data_sim.load_data()
        if data_sim.verify_name_in_data_names('label') == True:
            pass
        else:
            data_sim.insert_or_update_data({'label': np.array(['nuadm'+str(i)])})
        nuadm_sims.append(data_sim)
    
    plot_cum_oil(finescale_sim, nuadm_sims, 'cum_oil_ameba.svg')
    
    import pdb; pdb.set_trace()





    pass
