import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
import os

def run():


    fig, ax = plt.subplots(figsize=(6, 1), layout='constrained')

    # cmap = (mpl.colors.ListedColormap(['cyan', 'green', 'yellow'])
    #         .with_extremes(over='red', under='blue'))

    norm = mpl.colors.Normalize(vmin=1, vmax=13.14)

    # fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap='gist_rainbow_r'),
    #             cax=ax, orientation='horizontal', label='Some Units')

    colors = [(1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 1, 1), (0, 0, 1)]  # R -> G -> B
    colors.reverse()
    cmap_name = 'my_list'
    n_bin = 200
    cmap = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bin)
    ticks = [1, 13.14]
    
    cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
                cax=ax, orientation='horizontal', aspect=100)
    cbar.set_ticks(ticks)
    
    
    


    fig.savefig(os.path.join('results', 'test_colorbar.png'))