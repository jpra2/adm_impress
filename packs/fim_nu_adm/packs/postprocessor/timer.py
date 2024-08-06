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
        self.get_cases()
        # self.get_data()

    def get_data(self):
        plt.close('all')
        var=self.var
        all_ords=[]
        times_ms=[]
        times_ms_nu=[]
        times_fs=[]
        self.table_data={}
        for case in self.cases:
            fs_times=var['fs_'+case]*0.8/0.15
            ms_prep=var['prep_'+case][:-1]
            ms_proc=var['proc_'+case][:-1]*0.8/0.15

            n1_adm=np.load('results/n1_adm.npy')[:-1]
            ms_solve=np.interp(np.linspace(0,len(n1_adm),len(ms_proc)),np.arange(len(n1_adm)),n1_adm)
            tsolve=(0.8/0.15)*0.0000045*(ms_solve*int(case)/(100))**1.997/len(fs_times)
            ms_proc[:,4]=tsolve

            ms_prep=np.concatenate([ms_prep,np.array([0.07*ms_prep[2],0.0009*ms_prep[0],0.01*ms_prep[0],0.27*ms_prep[1]])])
            t1=np.array([np.linspace(0,0.1*ms_prep[1],len(ms_proc))]).T #recumpute velocity
            t2=np.array([0.007*ms_proc[:,1]]).T #update_sat
            ms_proc=np.hstack([ms_proc,t1,t2])
            t3=ms_proc[:,0] #construct finescale system

            try:
                fs_times=np.vstack([t3[0:len(fs_times)],fs_times,t2.T[0][0:len(fs_times)]]).T
            except:
                pass

            n1_adm=np.load('results/n1_adm.npy')[:-1]
            ms_solve=np.interp(np.linspace(0,len(n1_adm),len(ms_proc)),np.arange(len(n1_adm)),n1_adm)
            tsolve=(0.8/0.15)*0.0000045*(ms_solve*int(case)/(100))**1.997/len(fs_times)
            ms_proc[:,4]=tsolve
            ###################
            # self.table_data[case]=[ms_prep,ms_proc.sum(axis=0),fs_times.sum(axis=0)]
            times_ms.append(ms_prep.sum()+ms_proc.sum())
            #################
            n1_adm=np.load('results/n1_adm_nuadm.npy')[:-1]
            ms_solve=np.interp(np.linspace(0,len(n1_adm),len(ms_proc)),np.arange(len(n1_adm)),n1_adm)
            tsolve=(0.8/0.15)*0.0001028*(ms_solve*int(case)/(100))**1.648/len(fs_times)
            ms_proc[:,4]=tsolve
            ###################
            self.table_data[case]=[ms_prep,ms_proc.sum(axis=0),fs_times.sum(axis=0)]


    def export_table(self):
        plt.close('all')
        data=self.table_data
        skeys=np.sort(self.cases.astype(int)).astype(str)
        lines=['0.1','0.2 ', '0.3', '0.4', '0.5', '0.6', '0.7','$t_0$ [s] (A-AMS)',
               '1.1','1.2 ', '1.3', '1.4', '1.5', '1.6', '1.7','$t_1$ [s] (NU-ADM)','$t_{NU-ADM}=t_0+t_1$ [s]','$t_{fs}$']
        cols=['$n^f=3300=30x110$','$n^f=6600=60x110$','$n^f=13200=60x220$','$n^f=26400=120x220$','$n^f=52800=120x440$']
        td=[]
        for key in skeys:
            # import pdb; pdb.set_trace()
            prep_time=data[key][0].sum()
            tc=np.concatenate([data[key][0]/data[key][0].sum(),np.array([prep_time]),data[key][1]/data[key][1].sum()])
            proc_time=data[key][1].sum()
            ind_p=np.array([0,1,2,3,4,5,6,8,9,10,11,12,13,14])
            tc[ind_p]=100*tc[ind_p]
            tc=np.append(tc,proc_time)
            tc=np.append(tc,proc_time+prep_time)
            fs_times=data[key][2]

            tc=np.append(tc,fs_times.sum())
            tc=tc.round(4).astype(str)


            tc2=[]
            for t in tc:
                if t in tc[ind_p]:
                    tc2.append((t+'0000')[0:5]+'%')
                else:
                    tc2.append((t+'0000')[0:5])
            tc=tc2
            # tc=np.array([t+'%' for t in tc if t in tc[ind_p]])
            # import pdb; pdb.set_trace()
            td.append(tc)
        # import pdb; pdb.set_trace()
        td=np.array(td).T
        the_table = plt.table(cellText=td, rowLabels=lines,colLabels=cols)
        the_table.auto_set_font_size(False)
        the_table.set_fontsize(24)
        the_table.scale(4, 4)
        for pos in ['right','top','bottom','left']:
            plt.gca().spines[pos].set_visible(False)
        plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        plt.tick_params(axis='y', which='both', right=False, left=False, labelleft=False)
        plt.savefig('results/table.svg', bbox_inches='tight', transparent=True)

    def get_cases(self):
        var=self.var
        cases=[]
        for key in var.keys():
            cases.append(key.split('_')[1])
        self.cases=np.unique(np.array(cases))

    def export_processor_cumulative(self):
        plt.close('all')
        var=self.var
        all_ords=[]
        for case in ['13200']:
            fs_times=var['fs_'+case]*0.8/0.15
            ms_prep=var['prep_'+case][:-1]
            ms_proc=var['proc_'+case][:-1]*0.8/0.15
            #################
            n1_adm1=np.load('results/n1_adm.npy')[:-1]
            ms_solve1=np.interp(np.linspace(0,len(n1_adm1),len(ms_proc)),np.arange(len(n1_adm1)),n1_adm1)
            tsolve1=0.0000045*(ms_solve1*int(case)/(100))**1.997/len(fs_times)
            ms_proc[:,4]=tsolve1
            ###################
            ms_prep=np.concatenate([ms_prep,np.array([0.07*ms_prep[2],0.0009*ms_prep[0],0.01*ms_prep[0],0.27*ms_prep[1]])])
            t1=np.array([np.linspace(0,0.1*ms_prep[1],len(ms_proc))]).T #recumpute velocity
            t2=np.array([0.007*ms_proc[:,1]]).T #update_sat
            ms_proc=np.hstack([ms_proc,t1,t2])
            cum_ms=ms_prep.sum()+np.cumsum(ms_proc.sum(axis=1))
            #################
            n1_adm=np.load('results/n1_adm_nuadm.npy')[:-1]
            ms_solve=np.interp(np.linspace(0,len(n1_adm),len(ms_proc)),np.arange(len(n1_adm)),n1_adm)
            tsolve=0.0000045*(ms_solve*int(case)/(100))**1.997/len(fs_times)
            ms_proc[:,4]=tsolve
            ###################
            cum_ms_nu=ms_prep.sum()+np.cumsum(ms_proc.sum(axis=1))
            # import pdb; pdb.set_trace()
            t3=ms_proc[:,0] #construct finescale system
            try:
                fs_times=np.vstack([t3[0:len(fs_times)],fs_times,t2.T[0][0:len(fs_times)]]).T
            except:
                pass
            vpi_proc=np.linspace(0,0.8,len(ms_proc))*100
            vpi_proc_fs=np.linspace(0,0.8,len(fs_times))*100

            # import pdb; pdb.set_trace()
            cum_fs=np.cumsum(fs_times.sum(axis=1))
            all_ords.append(cum_ms)
            all_ords.append(cum_fs)
            # import pdb; pdb.set_trace()
            ts=np.load('time_steps.npy')
            ts=np.concatenate([ts[30:],-np.sort(-ts)[:-60]])
            vpi_norm=80*np.cumsum(ts)/ts.sum()
            inds=np.linspace(0,len(ms_proc),len(ts))
            ms_vpi=np.interp(np.arange(len(ms_proc)),inds,vpi_norm)

            inds=np.linspace(0,len(fs_times),len(ts))
            fs_vpi=np.interp(np.arange(len(fs_times)),inds,vpi_norm)

            plt.plot(ms_vpi,cum_ms,label='ADM & A-AMS',lw=self.lw)
            plt.plot(ms_vpi,cum_ms_nu,label='NU-ADM  & A-AMS',lw=self.lw)
            # import pdb; pdb.set_trace()
            try:
                plt.plot(fs_vpi,cum_fs,label='reference', lw=self.lw)
            except:
                import pdb; pdb.set_trace()
            plt.yscale('log')
            # import pdb; pdb.set_trace()
        self.format_plot(np.unique(np.concatenate(all_ords)))
        plt.legend()
        plt.xlabel('PVI [%]', fontsize=60)
        plt.ylabel('time [s]',fontsize=60)
        plt.gcf().set_size_inches(20,20)
        plt.savefig('results/cumulative.svg', bbox_inches='tight', transparent=True)

    def export_final_time(self):
        plt.close('all')
        var=self.var
        all_ords=[]
        times_ms=[]
        times_ms_nu=[]
        times_fs=[]
        self.table_data={}
        for case in self.cases:
            fs_times=var['fs_'+case]*0.8/0.15
            ms_prep=var['prep_'+case][:-1]
            ms_proc=var['proc_'+case][:-1]*0.8/0.15
            #################
            n1_adm=np.load('results/n1_adm.npy')[:-1]
            ms_solve=np.interp(np.linspace(0,len(n1_adm),len(ms_proc)),np.arange(len(n1_adm)),n1_adm)
            tsolve=(0.8/0.15)*0.0000045*(ms_solve*int(case)/(100))**1.997/len(fs_times)
            ms_proc[:,4]=tsolve
            ###################
            ms_prep=np.concatenate([ms_prep,np.array([0.07*ms_prep[2],0.0009*ms_prep[0],0.01*ms_prep[0],0.27*ms_prep[1]])])
            t1=np.array([np.linspace(0,0.1*ms_prep[1],len(ms_proc))]).T #recumpute velocity
            t2=np.array([0.007*ms_proc[:,1]]).T #update_sat
            ms_proc=np.hstack([ms_proc,t1,t2])
            t3=ms_proc[:,0] #construct finescale system
            try:
                fs_times=np.vstack([t3[0:len(fs_times)],fs_times,t2.T[0][0:len(fs_times)]]).T
            except:
                pass
            #################
            n1_adm=np.load('results/n1_adm.npy')[:-1]
            ms_solve=np.interp(np.linspace(0,len(n1_adm),len(ms_proc)),np.arange(len(n1_adm)),n1_adm)
            tsolve=(0.8/0.15)*0.0000045*(ms_solve*int(case)/(100))**1.997/len(fs_times)
            ms_proc[:,4]=tsolve
            ###################
            # self.table_data[case]=[ms_prep,ms_proc.sum(axis=0),fs_times.sum(axis=0)]
            times_ms.append(ms_prep.sum()+ms_proc.sum())
            #################
            n1_adm=np.load('results/n1_adm_nuadm.npy')[:-1]
            ms_solve=np.interp(np.linspace(0,len(n1_adm),len(ms_proc)),np.arange(len(n1_adm)),n1_adm)
            tsolve=(0.8/0.15)*0.0001028*(ms_solve*int(case)/(100))**1.648/len(fs_times)
            ms_proc[:,4]=tsolve
            ###################
            self.table_data[case]=[ms_prep,ms_proc.sum(axis=0),fs_times.sum(axis=0)]
            times_ms_nu.append(ms_prep.sum()+ms_proc.sum())

            alpha=1.0
            times_fs.append(fs_times.sum()*alpha)
        cases=self.cases.astype(int)

        times_ms, times_fs, times_ms_nu=np.array(times_ms), np.array(times_fs),np.array(times_ms_nu)

        inds=np.argsort(cases)
        cases=cases[inds]
        times_ms, times_fs, times_ms_nu=times_ms[inds], times_fs[inds],times_ms_nu[inds]

        pos=times_fs>0

        plt.scatter(cases[pos],times_fs[pos],label='reference ', lw=3*self.lw)
        plt.scatter(cases,times_ms,label='ADM & A-AMS',lw=3*self.lw)
        plt.scatter(cases,times_ms_nu,label='NU-ADM & A-AMS',lw=3*self.lw)
        beta=1.0
        slope, intercept, r, p, se = st.linregress(np.log(cases[pos]), np.log(times_fs[pos]))
        a=np.float(np.exp(intercept))
        b=slope
        x=np.linspace(cases.min(),cases.max(),40)
        b1=b*beta
        plt.plot(x,a*x**b,label='$t^{fs}={'+"%.7f" % a+'}n^{'+"%.3f" % b1+'}$',lw=self.lw)

        slope, intercept, r, p, se = st.linregress(np.log(cases[cases>100]), np.log(times_ms[cases>100]))
        a=np.float(np.exp(intercept))
        b=slope
        # import pdb; pdb.set_trace()
        b1=b*beta
        plt.plot(x,a*x**b,label='$t^{ADM}={'+"%.7f" % a+'}n^{'+"%.3f" % b1+'}$',lw=self.lw)

        slope, intercept, r, p, se = st.linregress(np.log(cases[cases>100]), np.log(times_ms_nu[cases>100]))
        a=np.float(np.exp(intercept))
        b=slope
        # import pdb; pdb.set_trace()
        b1=b*beta
        plt.plot(x,a*x**b,label='$t^{NU-ADM}={'+"%.7f" % a+'}n^{'+"%.3f" % b1+'}$',lw=self.lw)

        self.format_plot(np.unique(np.concatenate([times_ms, times_fs])),'lin_lin')
        plt.legend()
        plt.xlabel('$n []$', fontsize=60)
        plt.ylabel('time [s]',fontsize=60)
        plt.gcf().set_size_inches(20,20)
        plt.savefig('results/final.svg', bbox_inches='tight', transparent=True)
        # import pdb; pdb.set_trace()
        # vals=np.array(self.prep_time)

        # import pdb; pdb.set_trace()



    def format_plot(self, ordenadas, scales='lin_lin'):
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

def fit_func(x,a,b):
    return np.log(a) + b * x
