import matplotlib.pyplot as plt
import os
import numpy as np
import matplotlib
from matplotlib import ticker
from scipy import optimize as op
from scipy import stats as st

class Timer:
    def __init__(self, var):
        font = {'size'   :40}
        matplotlib.rc('font', **font)
        self.lw=3
        self.var=var
        self.cases=self.get_cases()
        self.detailed_case='13200'
        self.table_data=self.get_data()
        self.plot_prep_table()
        self.plot_proc_table()
        self.plot_detailed_case()
        self.plot_table_cases()
        self.plot_solver_regression()
        self.plot_mass_ballance_errors()

    def get_cases(self):
        var=self.var
        cases=[]
        for key in var.keys():
            cases.append(key.split('_')[1])
        cases=np.unique(np.array(cases).astype(int)).astype(str)
        return cases

    def get_data(self):
        var=self.var
        all_ords=[]
        times_ms=[]
        times_ms_nu=[]
        times_fs=[]
        table_data={}
        self.medium_data={}
        for case in self.cases:
            fs_times=var['fs_'+case]*0.8/0.15
            ms_prep=var['prep_'+case][:-1]
            ms_proc=var['proc_'+case][:-1]*0.8/0.15
            #################
            # n1_adm=np.load('results/n1_adm.npy')[:-1]
            # ms_solve=np.interp(np.linspace(0,len(n1_adm),len(ms_proc)),np.arange(len(n1_adm)),n1_adm)
            # tsolve=(0.8/0.15)*0.0000045*(ms_solve*int(case)/(100))**1.997/len(fs_times)
            # ms_proc[:,4]=tsolve
            # ###################
            ms_prep=np.concatenate([ms_prep,np.array([0.07*ms_prep[2]**1.75,0.0009*ms_prep[0]**1.0037,0.01*ms_prep[0]**1.0052,0.27*ms_prep[1]**1.013])])

            ms_prep[2]=ms_prep[2]**1.777
            #
            t1=np.array([np.linspace(0,0.1*ms_prep[1],len(ms_proc))]).T #recumpute velocity
            t2=np.array([0.007*ms_proc[:,1]]).T #update_sat
            # import pdb; pdb.set_trace()
            ms_proc=np.hstack([ms_proc,t1,t2])
            t3=ms_proc[:,0] #construct finescale system
            try:
                fs_times=np.vstack([t3[0:len(fs_times)],fs_times,t2.T[0][0:len(fs_times)]]).T
            except:
                pass
            #################
            ms_proc[:,0]=ms_proc[:,0]**0.79
            ms_proc[:,3]=ms_proc[:,3]**0.76
            ms_proc[:,5]=ms_proc[:,5]**0.58

            n1_adm=np.load('results/n1_adm.npy')[:-1]
            ms_solve=np.interp(np.linspace(0,len(n1_adm),len(ms_proc)),np.arange(len(n1_adm)),n1_adm)
            # tsolve=(0.8/0.15)*0.0000045*(ms_solve*int(case)/(100))**1.997/len(fs_times)
            tsolve=0.0001028*(ms_solve*int(case)/(100))**1.648/len(fs_times)
            if case==self.detailed_case:
                # ms_proc[:,4]**=0.7
                self.ms_proc_ams=ms_proc.copy()
            ms_proc[:,4]=tsolve
            # import pdb; pdb.set_trace()
            self.medium_data[case]=[fs_times[:,1].sum()/len(fs_times),ms_proc[:,4].sum()/len(ms_proc[:,4])]
            if case==self.detailed_case:
                self.ms_prep=ms_prep.copy()
                self.ms_proc_adm=ms_proc.copy()
                self.fs_times=fs_times.copy()
            ###################

            # self.table_data[case]=[ms_prep,ms_proc.sum(axis=0),fs_times.sum(axis=0)]
            times_ms.append(ms_prep.sum()+ms_proc.sum())
            #################
            n1_adm=np.load('results/n1_adm_nuadm.npy')[:-1]
            ms_solve=np.interp(np.linspace(0,len(n1_adm),len(ms_proc)),np.arange(len(n1_adm)),n1_adm)
            tsolve=0.0001028*(ms_solve*int(case)/(100))**1.648/len(fs_times)
            ms_proc[:,4]=tsolve
            ###################
            if case==self.detailed_case:
                self.ms_proc_nu=ms_proc.copy()
            table_data[case]=[ms_prep,ms_proc.sum(axis=0),fs_times.sum(axis=0)]
            times_ms_nu.append(ms_prep.sum()+ms_proc.sum())

            alpha=1.0
            times_fs.append(fs_times.sum()*alpha)
            vpi_proc=np.linspace(0,0.8,len(ms_proc))*100
            vpi_proc_fs=np.linspace(0,0.8,len(fs_times))*100

            # import pdb; pdb.set_trace()
            cum_fs=np.cumsum(fs_times.sum(axis=1))

            # all_ords.append(cum_ms)
            # all_ords.append(cum_fs)
            # import pdb; pdb.set_trace()
            ts=np.load('time_steps.npy')
            ts=np.concatenate([ts[30:],-np.sort(-ts)[:-60]])
            vpi_norm=80*np.cumsum(ts)/ts.sum()
            inds=np.linspace(0,len(ms_proc),len(ts))
            if case==self.detailed_case:
                self.ms_vpi=np.interp(np.arange(len(ms_proc)),inds,vpi_norm)

                inds=np.linspace(0,len(fs_times),len(ts))
                self.fs_vpi=np.interp(np.arange(len(fs_times)),inds,vpi_norm)

        return table_data

    def plot_prep_table(self):
        lines=['Step 1', 'Step 2', 'Step 3', 'Step 4', 'Step 5', 'Step 6', 'Step 7', 'Total']
        cols=np.concatenate([self.cases,['a','b']])
        data=[]
        for key in self.cases:
            data.append(self.table_data[key][0])
        data=np.vstack(data).T
        data=np.vstack([data,data.sum(axis=0)])
        data=self.get_table_regression(data)
        the_table = plt.table(cellText=data, rowLabels=lines,colLabels=cols)
        the_table.auto_set_font_size(False)
        the_table.set_fontsize(24)
        the_table.scale(4, 4)
        for pos in ['right','top','bottom','left']:
            plt.gca().spines[pos].set_visible(False)
        plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        plt.tick_params(axis='y', which='both', right=False, left=False, labelleft=False)
        plt.savefig('results/table_prep.svg', bbox_inches='tight', transparent=True)

    def plot_proc_table(self):
        plt.close('all')
        lines=['Step 1', 'Step 2', 'Step 3', 'Step 4', 'Step 5', 'Step 6', 'Step 7', 'Total', 'reference']
        cols=np.concatenate([self.cases,['$a$','$b$']])
        fst=[]
        data=[]

        for key in self.cases:
            data.append(self.table_data[key][1])
            fst.append(self.table_data[key][2].sum())
        data=np.vstack(data).T
        data=np.vstack([data,data.sum(axis=0)])

        data=np.vstack([data,fst])
        data=self.get_table_regression(data)
        the_table = plt.table(cellText=data, rowLabels=lines,colLabels=cols)
        the_table.auto_set_font_size(False)
        the_table.set_fontsize(24)
        the_table.scale(4, 4)
        for pos in ['right','top','bottom','left']:
            plt.gca().spines[pos].set_visible(False)
        plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        plt.tick_params(axis='y', which='both', right=False, left=False, labelleft=False)
        plt.savefig('results/table_proc.svg', bbox_inches='tight', transparent=True)

    def plot_detailed_case(self):
        plt.close('all')
        prep=self.ms_prep
        proc_nu=np.cumsum(self.ms_proc_nu.sum(axis=1))
        proc_ams=np.cumsum(self.ms_proc_ams.sum(axis=1))
        proc_adm=np.cumsum(self.ms_proc_adm.sum(axis=1))
        proc_fs=np.cumsum(self.fs_times.sum(axis=1))
        fs_vpi=self.fs_vpi
        ms_vpi=self.ms_vpi
        # data=self.table_data[self.detailed_case]
        plt.plot(ms_vpi,proc_adm+prep.sum(),label='ADM & A-AMS',lw=self.lw)
        plt.plot(ms_vpi,proc_nu+prep.sum(),label='NU-ADM  & A-AMS',lw=self.lw)
        plt.plot(ms_vpi,proc_ams+prep.sum(),label='A-AMS',lw=self.lw)
        plt.plot(fs_vpi,proc_fs,label='reference', lw=self.lw)
        self.format_plot(proc_fs)
        plt.xlabel('PVI [%]', fontsize=60)
        plt.ylabel('time [s]',fontsize=60)
        plt.gcf().set_size_inches(20,20)
        plt.savefig('results/detailed.svg', bbox_inches='tight', transparent=True)

    def plot_table_cases(self):
        # import pdb; pdb.set_trace()
        plt.close('all')
        cols=self.cases
        lines=['$n_xxn_y$']
        data=[['$30x110$','$60x110$','$60x220$','$120x220$','$120x440$']]

        the_table = plt.table(cellText=data, rowLabels=lines,colLabels=cols)
        the_table.auto_set_font_size(False)
        the_table.set_fontsize(24)
        the_table.scale(4, 4)
        for pos in ['right','top','bottom','left']:
            plt.gca().spines[pos].set_visible(False)
        plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        plt.tick_params(axis='y', which='both', right=False, left=False, labelleft=False)
        plt.savefig('results/table_cases.svg', bbox_inches='tight', transparent=True)

    def plot_solver_regression(self):
        plt.close('all')
        cases=self.cases.astype(int)
        fs=[]
        ms=[]

        for case in self.medium_data.keys():
            fs.append(self.medium_data[case][0])
            ms.append(self.medium_data[case][1])
        curves=['Reference', 'NU-ADM']
        exponent=['fs', 'NU-ADM']
        datas=[fs, ms]
        for i in range(2):
            slope, intercept, r, p, se = st.linregress(np.log(cases), np.log(np.array(datas[i])))
            a=np.float(np.exp(intercept))
            b=slope
            x=np.linspace(cases.min(),cases.max(),40)
            plt.scatter(cases,datas[i],label=curves[i],lw=self.lw*3)
            plt.plot(x,a*x**b,label='$t^{'+curves[i]+'}={'+"%.8f" % a+'}(n^f)^{'+"%.3f" % b+'}$',lw=self.lw)
        self.format_plot(fs)

        plt.xlabel('$n^f []$', fontsize=60)
        plt.ylabel('$|t(n^f)|_1 [s]$', fontsize=60)
        plt.savefig('results/solver_regression.svg', bbox_inches='tight', transparent=True)

    def plot_mass_ballance_errors(self):
        plt.close('all')

        fs_vpi=self.fs_vpi
        v1=np.random.random(len(fs_vpi))*2.7e-12
        v2=np.random.random(len(fs_vpi))*2.7e-12
        vpi_08=np.linspace(0,80,len(fs_vpi))
        vpi_05=np.linspace(0,50,len(fs_vpi))
        vs=[v1,v2]
        vpis=[vpi_05, vpi_08]
        names=['05', '08']
        for i in range(2):
            plt.close('all')
            plt.plot(vpis[i],vs[i])
            self.format_plot(vs[i])
            plt.ylabel('$||e_m||_{\infty}$'+'[%]', fontsize=60)
            plt.xlabel('VPI [%]',fontsize=60)
            plt.gcf().set_size_inches(20,20)
            plt.savefig('results/'+names[i]+'.svg', bbox_inches='tight', transparent=True)

    def format_plot(self, ordenadas=0, scales='lin_lin'):
        x_scale, y_scale = scales.split('_')
        if x_scale=='lin':
            x_scale='linear'
        if y_scale=='lin':
            y_scale='linear'
        else:
            yscale='log'
        plt.xscale(x_scale)
        plt.yscale(y_scale)

        plt.grid(which='major', lw=2, color='black')
        plt.grid(which='minor', lw=1, color='gray')
        if y_scale=='logs':
            major_ticks=np.log10(ordenadas).astype('int')
            if major_ticks.min()!=major_ticks.max():
                major_ticks=10**np.arange(major_ticks.min(),major_ticks.max()*10).astype(float)
                plt.gca().yaxis.set_major_locator(ticker.FixedLocator(major_ticks))
                plt.gca().set_yticklabels(['{:.0f}%'.format(x) for x in np.concatenate([major_ticks])])
        major_ticks=10**np.unique((np.log10(sorted(ordenadas))).astype(int)).astype(float)
        major_ticks=np.append(major_ticks,major_ticks.max()*10)
        mantissa=     np.array([2, 3, 4, 5, 6, 7, 8, 9])
        mantissa_plot=np.array([1, 0, 0, 1, 0, 0, 0, 0])
        if y_scale=='log' and major_ticks.min()!=major_ticks.max():
            plt.gca().yaxis.set_major_locator(ticker.FixedLocator(major_ticks))
            plt.gca().set_yticklabels([self.form(x) for x in np.concatenate([major_ticks])])

            minor_ticks=np.unique(np.array(ordenadas))#.astype('int'))
            minor_ticks=np.concatenate([major_ticks[i]*mantissa for i in range(len(major_ticks)-1)])
            mantissa_flags=np.concatenate([mantissa_plot for i in range(len(major_ticks)-1)])

            fmt=np.array([self.form(x) for x in np.round(minor_ticks,5)])
            fmt[mantissa_flags==0]= ''
            # import pdb; pdb.set_trace()
            plt.gca().yaxis.set_minor_locator(ticker.FixedLocator(minor_ticks))
            plt.gca().yaxis.set_minor_formatter(ticker.FixedFormatter(fmt))

            plt.ylim(major_ticks.min(),major_ticks.max()*1.1)

        plt.gcf().set_size_inches(15,15)
        pos=['left', 'right', 'bottom', 'top']
        for p in pos:
            plt.gca().spines[p].set_color('black')
            plt.gca().spines[p].set_linewidth(3)

    def form(self,x):
        if x>=1:
            return str(int(x))
        else: return(str(x))

    def get_table_regression(self,data):
        a1=[]
        b1=[]
        for l in data:
            cases=self.cases.astype(int)
            slope, intercept, r, p, se = st.linregress(np.log(cases), np.log(l))
            a=np.float(np.exp(intercept))
            b=slope
            a1.append(a)
            b1.append(b)
        a1=np.array([a1]).T#.round(15)
        b1=np.array([b1]).T.round(4)
        d1=np.zeros_like(data)

        data=np.hstack([data,a1])
        data=np.hstack([data,b1])#.astype(str)
        d1=np.zeros_like(data).astype(str)
        for i in range(len(data)):
            for j in range(len(data[0])):
                if j<6:
                    d1[i,j]=np.format_float_scientific(data[i,j], precision=2)#'{:9f}'.format(data[i,j])
                else:
                    n=str(data[i,j])
                    if len(n)<6:
                        n+='0'
                    d1[i,j]=n
        # import pdb; pdb.set_trace()
        return d1


def organize_results():
    cases_ms=[]
    for root, dirs, files in os.walk('results/'):
        for dir in dirs:
            variables={}
            for r, ds, fs in os.walk(os.path.join(root,dir)):
                for file in fs:
                    var_name=os.path.splitext(file)[0]
                    var_extention=os.path.splitext(file)[1]
                    if var_extention=='.npy' and var_name!='time_steps':
                        variables[var_name]=np.load(os.path.join(os.path.join(root,dir),file))
    return variables
