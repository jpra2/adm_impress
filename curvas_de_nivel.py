import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt

# delta = 1.0
# x = np.arange(-3.0, 3.0, delta)
# y = np.arange(-3.0, 3.0, delta)
# x, y = np.meshgrid(x, y)
#
# z = np.sqrt(9-x**2+7*y**2)



def plot_curva(x,y,z, number=5):
    plt.close('all')
    ax = plt.axes(projection='3d')
    #zline = np.linspace(0, 15, 1000)
    #xline = np.sin(zline)
    #yline = np.cos(zline)
    ax.scatter3D(x, y, np.log10(z), 'gray')

    #
    # ind=np.argsort(x)
    # nx=(x[ind]==x[ind[0]]).sum()
    # ny=int(len(ind)/nx)
    # x=x[ind].reshape(ny,nx)
    # y=y[ind].reshape(ny,nx)
    # z=z[ind].reshape(ny,nx)
    #
    # fig, ax = plt.subplots()
    # CS = ax.contour(x, y, z,number,colors='k')
    # ax.clabel(CS, inline=0.5, fontsize=8)
    # ax.set_title('Curvas de nivel')
    plt.savefig('results/curve_img.png',bbox_inches='tight')
    plt.tight_layout()
    import pdb; pdb.set_trace()
# plot_curva(x,y,z, number=10)
