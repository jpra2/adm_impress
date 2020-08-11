import numpy as np
import os
import matplotlib.pyplot as plt
import shutil
import matplotlib
from matplotlib import ticker
from brokenaxes import brokenaxes


font = {'size'   : 40}
matplotlib.rc('font', **font)
def remove_previous_files():
    # Cleans the previous results
    shutil.rmtree('results/biphasic/ms')
    shutil.rmtree('results/biphasic/finescale')
    # creates new folders
    os.mkdir('results/biphasic/ms')
    os.mkdir('results/biphasic/finescale')
    os.mkdir('results/biphasic/finescale/vtks')

def run_test_cases():
    vpis_for_save=np.array([0.0])

    np.save('flying/vpis_for_save.npy',vpis_for_save)
    # os.system("python testting2_biphasic.py")
    crs=[[3,7,1],[5, 11, 1],[9, 19, 1],[13, 27, 1]]
    np.save('flying/all_crs.npy',np.array(crs))
    # neta_lim_dual_values=     [ np.inf,    0.5,    1.0,    2.0,   10.0,  100.0]#,    2.0,    5.0,   10.0,  100.0,   500.0,1000.0, np.inf]
    neta_lim_dual_values=     [    0.5,    1.0,   2.0,   10.0,   100.0, np.inf]#,  np.inf,    1.0]# np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf]#, np.inf, np.inf, np.inf]#,    2.0,    5.0,   10.0,  100.0,   500.0,1000.0, np.inf]
    neta_lim_finescale_values=[ np.inf, np.inf, np.inf, np.inf, np.inf, np.inf]#,  np.inf, np.inf]# np.inf,    0.5,    1.0,    5.0,    1.0,   10.0,    0.5,    1.0,    5.0,    1.0,   10.0,    0.5,    1.0,    5.0,    1.0,   10.0]#, np.inf, np.inf, np.inf]#, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf]
    type_of_refinement_values=[  'uni',  'uni',  'uni',  'uni',  'uni',  'uni']#,    'uni', 'uni']#  'uni',  'uni',  'uni',  'uni',  'uni',  'uni',  'uni',  'uni',  'uni',  'uni',  'uni',  'uni',  'uni',  'uni',  'uni',  'uni']#,  'uni',  'uni',  'uni']#,  'uni',  'uni',  'uni',  'uni',  'uni',  'uni',  'uni']
    phiK_raz_lim_values=      [ np.inf, np.inf, np.inf, np.inf, np.inf, np.inf]#,  np.inf, np.inf]# np.inf, np.inf, np.inf, np.inf, np.inf, np.inf,    3.0,    3.0,    3.0,    3.0,    3.0,   10.0,   10.0,   10.0,   10.0,   10.0]#, np.inf, np.inf, np.inf]#, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf]
    delta_sat_max=            [    1.1,    1.1,    1.1,    1.1,    1.1,    1.1]#,     1.1,    1.1]#    1.1,    1.1,    1.1,    1.1,    1.1,    1.1,    1.1,    1.1,    1.1,    1.1,    1.1,    1.1,    1.1,    1.1,    1.1,    1.1]#,    1.1,    1.1,    1.1]#,    1.1,    1.1,    1.1,    1.1,    1.0,    1.0,    1.0]
    cr_inds=                  [      3,      3,      3,      3,      3,      3]#,       1,      1]
    for i in range(len(neta_lim_dual_values)):
        np.save('flying/delta_sat_max.npy',np.array([delta_sat_max[i]]))
        np.save('flying/neta_lim_finescale.npy',np.array([neta_lim_finescale_values[i]]))
        np.save('flying/neta_lim_dual.npy',np.array([neta_lim_dual_values[i]]))
        np.save('flying/type_of_refinement.npy',np.array([type_of_refinement_values[i]]))
        np.save('flying/phiK_raz_lim.npy',np.array([phiK_raz_lim_values[i]]))
        np.save('flying/crs.npy',np.array(crs[cr_inds[i]]))
        ms_case='neta_'+str(neta_lim_dual_values[i])+'_alpha_'+str(neta_lim_finescale_values[i])+\
        '_type_'+type_of_refinement_values[i]+'_beta_'+str(phiK_raz_lim_values[i])+\
        '_delta_'+str(delta_sat_max[i])+'_CR_'+str(cr_inds[i])+'/'
        os.makedirs('results/biphasic/ms/'+ms_case,exist_ok=True)
        os.makedirs('results/biphasic/ms/'+ms_case+'vtks',exist_ok=True)
        np.save('flying/ms_case.npy',np.array([ms_case]))
        os.system("python ADM_b.py")

def organize_results():
    cases_ms=[]
    for root, dirs, files in os.walk("results/biphasic/ms"):
        for dir in dirs:
            variables={}

            for r, ds, fs in os.walk(os.path.join(root,dir)):
                for file in fs:
                    var_name=os.path.splitext(file)[0]
                    var_extention=os.path.splitext(file)[1]
                    if var_extention=='.npy':
                        variables[var_name]=np.load(os.path.join(os.path.join(root,dir),file))
            if dir!='vtks':
                cases_ms.append([dir,variables])
    finescale={}
    for r, ds, fs in os.walk("results/biphasic/finescale"): #print(fs)
        for file in fs:
            var_name=os.path.splitext(file)[0]
            var_extention=os.path.splitext(file)[1]
            if var_extention=='.npy':
                finescale[var_name]=np.load(os.path.join(r,file))

    case_finescale=[['finescale',finescale]]
    all_cases_results=np.array(cases_ms+case_finescale)

    return all_cases_results

def get_axis_lims(ordenada,max_raz,pre_inds=[]):
    ind_superates=[]
    if (ordenada==0).sum()==0:
        # reasons=ordenada[:-1]/ordenada[1:]
        # reasons[1/reasons>reasons]=1/reasons[1/reasons>reasons]
        reasons=abs(ordenada[:-1]-ordenada[1:])

        reasons=np.concatenate([[1],reasons])
        ind_superates=np.arange(len(ordenada))[(reasons>max_raz) | (reasons<1/max_raz)]
        ind_separates=np.unique(np.concatenate([pre_inds, ind_superates]))
        i_prev=0
        ylims=[]

        for ind in ind_superates:
            vals=np.array([ordenada[i_prev],ordenada[ind-1]])

            if vals.max()<1.1:
                ylims.append((vals.min()-0.006,vals.max()+0.006))
            elif vals.min()%1==0:
                ylims.append((int(vals.min())-0.2,int(vals.max())+0.2))
            else:
                ylims.append((int(vals.min())-1,int(vals.max())+1.01))
            i_prev=ind
        # if i_prev<=len(ordenada)-1:
        if ordenada[-1]<1.1:
            ylims.append((ordenada[-1]-0.008,ordenada[-1]+0.008))
        elif ordenada.max()%1==0:
            ylims.append((int(ordenada[-1])-0.2,int(ordenada[-1])+0.2))
        elif len(ind_superates)>0:
            ylims.append((int(ordenada[-1])-1,int(ordenada[-1])+1.01))
        else:
            ylims=[(ordenada.min(),ordenada.max())]
    else:
        if ordenada.min()!=ordenada.max():
            ylims=[(ordenada.min(),ordenada.max())]
        else:
            ylims=[(0,1)]

    return ylims, ind_superates



def merge_intervals(intervals):
    s = sorted(intervals, key=lambda t: t[0])
    m = 0
    for  t in s:
        if t[0] > s[m][1]:
            m += 1
            s[m] = t
        else:
            s[m] = (s[m][0], t[1])
    return s[:m+1]


def setup(ax, title='teste'):
    """Set up common parameters for the Axes in the example."""
    # only show the bottom spine
    ax.yaxis.set_major_locator(ticker.NullLocator())
    ax.spines['right'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['top'].set_color('none')

    # define tick positions
    ax.xaxis.set_major_locator(ticker.MultipleLocator(1.00))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.25))

    ax.xaxis.set_ticks_position('bottom')
    ax.tick_params(which='major', width=1.00, length=5)
    ax.tick_params(which='minor', width=0.75, length=2.5, labelsize=10)
    ax.set_xlim(0, 5)
    ax.set_ylim(0, 1)
    ax.text(0.0, 0.2, title, transform=ax.transAxes,
            fontsize=14, fontname='Monospace', color='tab:blue')

def print_results(all_cases):
    units={'vpi':'vpi [%]','wor':'wor []','t_comp':'comp_time [s]','delta_t':'time-step []',
        'n1_adm':'Nadm/Nf[%]','el2':r'$||e_p||_2$'+' [%]','elinf':r'$||e_p||_\infty$'+' [%]', 'es_L2':'es_L2 [%]',
        'es_Linf':'es_Linf [%]', 'vpis_for_save':'vpi [%]','coupl':"Percentage of enrichment [%]",
        'refinement':"Percentage at fine-scale [%]", 'ep_haji_L2':'ep_rel_ms_L2[]',
        'ep_haji_Linf':'ep_rel_ms_Linf[]','er_L2':r'$||e_r||_2$'+' []','er_Linf':r'$||e_r||_\infty$'+' []',
        'ev_L2':r'$||e_v||_2$'+' [%]','ev_Linf':r'$||e_v||_\infty$'+' [%]'}
    variables=all_cases[0][1].keys()
    single_vars={}
    for u in units.keys():
        single_vars[u]=[]
    single_vars['neta']=[]
    single_vars['alpha']=[]
    single_vars['beta']=[]
    single_vars['delta']=[]
    single_vars['CR']=[]
    names_single_vars=['neta', 'beta', 'alpha', 'delta', 'CR']

    for variable in variables:
        plt.close('all')
        ymin=np.inf
        ymax=-np.inf

        for case in all_cases:
            case_name=case[0]
            case_data=case[1]

            if case_name!='finescale':
                single_vars[variable].append(case_data[variable][0])
                # case name defined variables
                if variable=='el2':
                    # import pdb; pdb.set_trace()
                    input_params=case_name.split('_')
                    for i, param in enumerate(input_params):
                        if param in names_single_vars:

                            if input_params[i+1]=='inf':
                                single_vars[param].append(1000.0)
                            else:
                                single_vars[param].append(float(input_params[i+1]))


                # if case_name[5:8]!='inf':
                #     if variable=='el2':
                #         single_vars['neta'].append(float(case_name[5:8]))
                # else:
                #     if variable=='el2':
                #         single_vars['neta'].append(1000.0)

            # case_data['vpi']
            if case_name[5]!='i':
                style='-.'
            else:
                style='-'
            if variable!="vpi" and variable!='vpis_for_save':

                if case_name=='finescale':
                    if variable!="n1_adm" and variable!="el2" and variable!="elinf" and variable!='es_L2'and variable!='es_Linf' and variable!='refinement' and variable!='coupl':

                        try:
                            plt.plot(100*case_data['vpi'], case_data[variable],label=case_name)
                        except:
                            pass
                else:
                    if variable=='n1_adm':
                        plt.plot(100*case_data['vpi'], 100*case_data[variable][:-1]/case_data[variable][-1], style, label=case_name)
                        plt.xlabel(units['vpi'])
                        plt.ylabel(units[variable])
                    elif variable!='es_L2' and variable!='es_Linf':
                        if case_data[variable].min()<ymin:
                            ymin=case_data[variable].min()
                        if case_data[variable].max()>ymax:
                            ymax=case_data[variable].max()

                        plt.plot(100*case_data['vpi'], case_data[variable],style,label=case_name)

                        plt.xlabel(units['vpi'])
                        plt.ylabel(units[variable])

                    else: # Saturation variables
                        if case_data[variable].min()<ymin:
                            ymin=case_data[variable].min()
                        if case_data[variable].max()>ymax:
                            ymax=case_data[variable].max()

                        plt.plot(100*case_data['vpis_for_save'], case_data[variable],style,label=case_name)

                        plt.xlabel(units['vpis_for_save'])
                        plt.ylabel(units[variable])

        if variable!='vpi' and variable!='vpis_for_save':
            if (variable=='el2' or variable=='elinf' or variable=='es_L2' or variable=='es_Linf') and False:
                if ymin>0:
                    ymin=10**int(np.log10(ymin))
                    if ymin<1e-2:
                        ymin=1e-2
                if ymax>0:
                    ymax=10**int(np.log10(ymax)+1)
                # plt.yscale('log')
                y_ticks=plt.gca().get_yticks()[:-1]
                y_ticks=np.concatenate(np.vstack([y_ticks,np.append(y_ticks[1:]/2,y_ticks[-1])]).T)[:-1]
                plt.gca().set_yticks(y_ticks)
                ticks=plt.gca().get_yticks()

                pos_ymax=(ticks<ymax).sum()
                t_ymin=ticks[ticks>=ymin].min()

                t_ymax=ticks[pos_ymax]
                plt.yticks(ticks)

                plt.gca().set_yticklabels(['{:.0f}%'.format(x) for x in plt.gca().get_yticks()])
                plt.ylim(t_ymin,t_ymax)

            plt.gca().tick_params(which='minor', length=10)
            plt.gca().tick_params(which='major', length=15)

            plt.grid()
            plt.legend()
            plt.gcf().set_size_inches(20,20)
            plt.savefig('results/biphasic/'+variable+'.png')

    for var in single_vars:
        abcissa_var='neta'
        control_variable='CR'

        all_abcissa=np.array(single_vars[abcissa_var])
        pos=all_abcissa<10000
        all_abcissa=all_abcissa[pos]
        all_ordenada=np.array(single_vars[var])[pos]
        control_parameters=single_vars[control_variable]
        # import pdb; pdb.set_trace()
        ordenadas=[]
        abcissas=[]
        control_vals=np.unique(control_parameters)

        for p in np.unique(control_parameters):
            orden=all_ordenada[control_parameters==p]
            abcis=all_abcissa[control_parameters==p]

            if abcis.max()<1000:
                orden=np.concatenate([[all_ordenada[all_abcissa==1000][0]],orden])
                abcis=np.concatenate([[1000],abcis])
            ordenadas.append(orden)
            abcissas.append(abcis)

        plt.close('all')

        linear_yaxis=['refinement','coupl', 'ev_L2','ev_Linf']
        if var in ['elinf','el2', 'ev_L2', 'ev_Linf', 'refinement', 'coupl']:
            plt.close('all')
            fig=plt.figure()
            plt.grid()
            plt.xlabel(r'$\epsilon$ []',labelpad=-20)
            if var not in ['neta', 'alpha', 'beta', 'delta']:
                plt.ylabel(units[var])#, labelpad=-20)

            control=0

            for abcissa, ordenada in zip(abcissas, ordenadas):
                plt.gca().tick_params(which='minor', length=10)
                plt.gca().tick_params(which='major', length=15)
                ind_sort=np.argsort(abcissa)
                abcissa=abcissa[ind_sort]
                ordenada=ordenada[ind_sort]

                # if ordenada.min()>0:
                plt.xscale('log')
                if var not in linear_yaxis:
                    plt.yscale('log')

                # f, (ax, ax2) = plt.subplots(2, 1, sharey=True)

                # ax.set_xlim(abcissa.min(), 10.0)  # outliers only
                # ax2.set_xlim(999, 1001)  # most of the data
                abc=abcissa[abcissa!=30]
                ord=ordenada[abcissa!=30]
                if control_variable=='beta':
                    label=r'$\beta^\lim = {}$'
                    control_val=control_vals[control]
                else:
                    label='CR = {}'
                    all_crs=np.load('flying/all_crs.npy')
                    control_val = int(control_vals[control])
                    control_val = str(all_crs[control_val][0])+', '+str(all_crs[control_val][1])
                if control in [1,2]:
                    plt.plot(abc, ord,label=label.format(control_val),marker='o',markersize=15, linewidth=5)
                # plt.scatter(abcissa, ordenada, s=300)

                plt.grid(axis='both', which='major', ls='-',lw=3)
                plt.grid(axis='both', which='minor', ls='--', alpha=0.4, lw=2)
                positions = np.unique(all_abcissa.astype(int))
                positions=positions[positions!=30]
                labels = positions.astype(str)


                plt.gca().xaxis.set_major_locator(ticker.FixedLocator(positions))
                plt.gca().xaxis.set_major_formatter(ticker.FixedFormatter(labels))

                if var!='refinement' and var!='coupl':
                    positions = 10**np.arange(int(np.log10(all_ordenada).min()), int(np.log10(all_ordenada).max())+1)

                else:
                    positions = 10**np.array([0.0, 1.0])

                pp=[]
                ppp=[]
                for p in positions:
                    for i in range(2,10):
                        pp.append(i*p)
                    for i in range(3,8,2):
                        ppp.append(i*p)

                # ppp.append(int(all_ordenada.min())+1)
                if var in linear_yaxis:
                    pp.append(int(all_ordenada.min()))
                    pp.append(int(all_ordenada.max()))


                close_values=[]
                # import pdb; pdb.set_trace()
                for n in np.sort(np.unique(all_ordenada)):
                    close_values.append(min(ppp, key=lambda x:abs(x-n)))
                all_ticks=np.zeros_like(pp)
                close_values= np.unique(np.array(close_values))
                for value in close_values:
                    ind=pp==value
                    all_ticks[ind]=np.array(pp)[ind][0]
                at=[]

                for v in all_ticks:
                    if v==0:
                        at.append('')
                    else:
                        at.append(int(v))


                all_ticks=at
                if var in linear_yaxis:
                    positions=np.array([])
                    labels=np.array([])
                    if var=='refinement':
                        pp.append(1)
                        pp.append(10)
                    if var=='coupl':
                        pp=np.arange(0,np.array(pp).max()+5,5)
                    all_ticks=np.array(pp).astype(int)                

                labels = positions.astype(int).astype(str)
                plt.gca().yaxis.set_major_locator(ticker.FixedLocator(positions))
                plt.gca().yaxis.set_major_formatter(ticker.FixedFormatter(labels))

                plt.gca().yaxis.set_minor_locator(ticker.FixedLocator(pp))
                plt.gca().yaxis.set_minor_formatter(ticker.FixedFormatter(all_ticks))
                plt.legend()
                plt.gca().xaxis.set_minor_formatter(ticker.NullFormatter())
                # ticks=plt.gca().get_yticks()
                control+=1
            # plt.xscale('log')

            plt.gcf().set_size_inches(15,15)
            plt.savefig('results/single_phase/'+var+'.pdf', bbox_inches='tight')

    import pdb; pdb.set_trace()
