import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
from math import pi
from packs.utils.utils_old import get_box
flying = 'flying'
name = 'results'
arquivos = os.listdir(flying)
for  arq in arquivos:
    if  arq.startswith(name):
        datas = np.load('flying/results_Dietz_30dg_case_1288.npy', allow_pickle=True)

        x425_all = np.zeros([12,10])
        b=0
        import pdb; pdb.set_trace()
        for data in datas[1:]:
            Sw = data[5]
            centroids = data[6]
            inds = np.zeros([10,100])
            Sw_org = np.zeros([10,100])
            x425 = np.zeros(10)
            for i in range(0,10):
                for j in range(0,100):
                    p0 = [0.3048*j, 0.0, -0.3048*(i+1) - 0.175976362*j]
                    p1 = [0.3048*(j+1), 3.048, -0.3048*i - 0.175976362*(j+1)]
                    inds[i,j] = get_box(centroids,np.array([p0,p1]))
                inds = inds.astype(int)
                Sw_org[i,:] = Sw[inds[i,:]]
                ind_to_ind = np.argwhere(np.round(Sw_org[i,:],3)>0.4).ravel()
                ind = inds[i,ind_to_ind][0:2]
                #import pdb; pdb.set_trace()
                if Sw[ind[0]]>0.425:
                    ind_to_ind = np.argwhere(np.round(Sw_org[i,:],3)<0.4).ravel()
                    ind[1] = ind[0]
                    ind[0] = inds[i,ind_to_ind][-1]

                Slim = Sw[ind]
                x = centroids[ind][:,0]
                x425[i] = x[0] + (0.425 - Slim[0])/(Slim[1] - Slim[0])*(x[1] - x[0])

            x425_all[b,:] = 100 - x425/0.3048#/np.cos(pi/6)
            z_plot = np.linspace(-0.5,-9.5,10)
            b+=1
        x_ans = np.array([[15.025, 17.0467, 19.0037, 20.5887, 22.0605, 23.2088, 24.9232, 26.0877, 27.1713, 29.0312],
                          [38.994, 40.1585, 41.6141, 42.9888, 44.1695, 45.6898, 46.9513, 48.5848, 50.0243, 51.0108],
                          [61.006, 62.8012, 63.9981, 65.0817, 66.7961, 68.0091, 69.772, 71.0335, 72.1172, 73.88],
                          [83.7943, 84.9749, 86.6731, 87.967, 89.455, 90.4739, 91.9942, 93.1748, 94.8892, 95.9567]])
        
        x425_plot = np.zeros([4,10])
        x425_plot[0,:] = x425_all[1,:]
        x425_plot[1,:] = x425_all[4,:]
        x425_plot[2,:] = x425_all[7,:]
        x425_plot[3,:] = x425_all[10,:]
        plt.figure(1)
        plt.plot(x425_plot[0,:],z_plot, 'r',label = ('vpi = 0.1'))
        plt.plot(x425_plot[1,:],z_plot,'g', label = ('vpi = 0.2'))
        plt.plot(x425_plot[2,:],z_plot, 'b', label = ('vpi = 0.3'))
        plt.plot(x425_plot[3,:],z_plot, 'y', label = ('vpi = 0.4'))#, 'r', x_ans, z_plot, 'b')
        plt.plot(x_ans[0,:],z_plot, 'k', x_ans[1,:],z_plot, 'k', x_ans[2,:],z_plot, 'k')
        plt.plot(x_ans[3,:],z_plot, 'k', label='UTCOMP')
        plt.grid()
        plt.title('Dietz -30° Solution Example', y=1.08)
        plt.legend(('vpi = 0.1', 'vpi = 0.2', 'vpi = 0.3', 'vpi = 0.4','Benchmark'))
        plt.legend(bbox_to_anchor=(.48, 1.08), loc=9, borderaxespad=0., ncol = 5)
        plt.ylabel('Distance in Z-Direction (ft)')
        plt.xlabel('Distance in X-Direction (ft)')
        plt.savefig('results/compositional/profile_Dietz_comparison.png')
        import pdb; pdb.set_trace()
