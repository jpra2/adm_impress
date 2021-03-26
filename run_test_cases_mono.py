import numpy as np
import os
import matplotlib.pyplot as plt
import shutil
import matplotlib
from matplotlib import ticker



font = {'size'   : 30}
matplotlib.rc('font', **font)
folder='biphasic'
# folder='single_phase_L1'
np.save('flying/folder.npy',np.array([folder]))

def remove_previous_files():
    # Cleans the previous results
    shutil.rmtree('results/'+folder+'/ms')
    shutil.rmtree('results/'+folder+'/finescale')
    # creates new folders
    os.mkdir('results/'+folder+'/ms')
    os.mkdir('results/'+folder+'/finescale')
    os.mkdir('results/'+folder+'/finescale/vtks')

def run_test_cases():
    # vpis_for_save=np.arange(0.0,1.001,0.01)
    vpis_for_save=np.array([0.0,0.00001])
    np.save('flying/vpis_for_save.npy',vpis_for_save)
    # os.system('python testting2_biphasic.py')

    crs=[[9,9,1],[5, 11, 1],[9, 19, 1],[13, 27, 1]]
    np.save('flying/all_crs.npy',np.array(crs))
    neta_lim_dual_values=     [    5.0]
    neta_lim_finescale_values=[ np.inf]
    type_of_refinement_values=[  'uni']
    phiK_raz_lim_values=      [ np.inf]
    delta_sat_max=            [ np.inf]
    cr_inds=                  [      0]
    for i in range(len(neta_lim_dual_values)):
    # if False:
        np.save('flying/delta_sat_max.npy',np.array([delta_sat_max[i]]))
        np.save('flying/neta_lim_finescale.npy',np.array([neta_lim_finescale_values[i]]))
        np.save('flying/neta_lim_dual.npy',np.array([neta_lim_dual_values[i]]))
        np.save('flying/type_of_refinement.npy',np.array([type_of_refinement_values[i]]))
        np.save('flying/phiK_raz_lim.npy',np.array([phiK_raz_lim_values[i]]))
        np.save('flying/crs.npy',np.array(crs[cr_inds[i]]))
        ms_case='neta_'+str(neta_lim_dual_values[i])+'_alpha_'+str(neta_lim_finescale_values[i])+\
        '_type_'+type_of_refinement_values[i]+'_beta_'+str(phiK_raz_lim_values[i])+\
        '_delta_'+str(delta_sat_max[i])+'_CR_'+str(cr_inds[i])+'/'
        os.makedirs('results/'+folder+'/ms/'+ms_case,exist_ok=True)
        os.makedirs('results/'+folder+'/ms/'+ms_case+'vtks',exist_ok=True)
        np.save('flying/ms_case.npy',np.array([ms_case]))
        os.system('python ADM_b.py')

def organize_results():
    cases_ms=[]
    for root, dirs, files in os.walk('results/'+folder+'/ms'):
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
    for r, ds, fs in os.walk('results/'+folder+'/finescale'): #print(fs)
        for file in fs:
            var_name=os.path.splitext(file)[0]
            var_extention=os.path.splitext(file)[1]
            if var_extention=='.npy':
                finescale[var_name]=np.load(os.path.join(r,file))

    case_finescale=[['finescale',finescale]]
    all_cases_results=np.array(cases_ms+case_finescale)

    return all_cases_results


def print_results(all_cases):
    units={'vpi':'PVI [%]','wor':'wor []','t_comp':'comp_time [s]','delta_t':'time-step []',
        'n1_adm':'Na-adm/Nf[%]','el2':r'$||e_p||_2$'+' [%]','elinf':r'$||e_p||_\infty$'+' [%]', 'es_L2':'es_L2 [%]',
        'es_Linf':'es_Linf [%]', 'vpis_for_save':'PVI [%]','coupl':'Percentage of enrichment [%]',
        'refinement':'Percentage at fine-scale [%]', 'ep_haji_L2':'ep_rel_ms_L2[]',
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
        #all_cases[:,0]
        # import pdb; pdb.set_trace()
        # for case in all_cases[np.array([0,1,2,3,4,5,8,9,10])]:
        for case in all_cases[np.array([1,4,-1])]:
        # import pdb; pdb.set_trace()
        # for case in all_cases:
            case_name=case[0]
            case_data=case[1]
            if case_name!='finescale':
                single_vars[variable].append(case_data[variable][0])
                if variable=='el2':
                    input_params=case_name.split('_')
                    for i, param in enumerate(input_params):
                        if param in names_single_vars:

                            if input_params[i+1]=='inf':
                                single_vars[param].append(1000.0)
                            else:
                                single_vars[param].append(float(input_params[i+1]))



            if case_name[5]!='i':
                style='-.'
            else:
                style='-'
            if variable!='vpi' and variable!='vpis_for_save':

                if case_name=='finescale':
                    if variable!='n1_adm' and variable!='el2' and variable!='elinf' and variable!='es_L2'and variable!='es_Linf' and variable!='refinement' and variable!='coupl':

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
                        try:
                            plt.plot(100*case_data['vpi'], case_data[variable],style,label=case_name)
                        except:
                            pass

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
            plt.savefig('results/'+folder+'/'+variable+'.svg', bbox_inches=12)

    for var in single_vars:
        abcissa_var='alpha'
        control_variable='delta'

        all_abcissa=np.array(single_vars[abcissa_var])
        pos=all_abcissa<10000
        all_abcissa=all_abcissa[pos]
        all_ordenada=np.array(single_vars[var])[pos]
        control_parameters=single_vars[control_variable]

        ordenadas=[]
        abcissas=[]
        control_vals=np.unique(control_parameters)

        for p in np.unique(control_parameters):
            orden=all_ordenada[control_parameters==p]
            abcis=all_abcissa[control_parameters==p]

            if abcis.max()<1000:
                if (all_abcissa==1000).sum()>0:
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
                try:
                    ordenada=ordenada[ind_sort]
                except:
                    ordenada=np.repeat(ordenada[0],len(ind_sort))

                # plt.xscale('log')
                # if var not in linear_yaxis:
                #     plt.yscale('log')

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

                plt.grid(axis='both', which='major', ls='-',lw=7)
                plt.grid(axis='both', which='minor', ls='--', alpha=0.4, lw=7)
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

                if var in linear_yaxis:
                    pp.append(int(all_ordenada.min()))
                    pp.append(int(all_ordenada.max()))


                close_values=[]
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
                        # positions.append(.5)
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
            plt.savefig('results/single_phase/'+var+'.svg', bbox_inches='tight')


def print_results_2(all_cases):

    units={'vpi':'PVI [%]','wor':'wor []','t_comp':'comp_time [s]','delta_t':'time-step []',
        'n1_adm':'Na-adm/Nf[%]','el2':r'$||e_p||_2$'+' [%]','elinf':r'$||e_p||_\infty$'+' [%]', 'es_L2':'es_L2 [%]',
        'es_Linf':'es_Linf [%]', 'vpis_for_save':'PVI [%]','coupl':'Percentage of enrichment [%]',
        'refinement':'Percentage at fine-scale [%]', 'ep_haji_L2':'ep_rel_ms_L2[]',
        'ep_haji_Linf':'ep_rel_ms_Linf[]','er_L2':r'$||e_r||_2$'+' []','er_Linf':r'$||e_r||_\infty$'+' []',
        'ev_L2':r'$||e_v||_2$'+' [%]','ev_Linf':r'$||e_v||_\infty$'+' [%]'}
    variables=all_cases[0][1].keys()
    single_vars={}
    for variable in variables:
        plt.close('all')
        ymin=np.inf
        ymax=-np.inf
        cases4print=['finescale',
                    'neta_inf_alpha_inf_n_uni_beta_inf_delta_0.1s',
                    'neta_1.0_alpha_1.0_uni_beta_3_delta_0.1',
                    'neta_1.0_alpha_1.0_uni_beta_3_delta_2.0',
                    'neta_1.0_alpha_1.0_uni_beta_inf_delta_0.1s',
                    'neta_inf_alpha_1.0_uni_beta_inf_delta_0.1s',
                    'neta_inf_alpha_1.0_uni_beta_3_delta_0.1s']


        names4print={cases4print[0]:'Reference',
                    cases4print[1]:'Traditional ADM',
                    # cases4print[2]:'A-ADM+Enrichment',
                    cases4print[2]:r'$ \delta = 0.1$',
                    # cases4print[3]:'Pure A-ADM+Enrichment',
                    cases4print[3]:r'$\delta = \infty$',
                    cases4print[4]:'Adaptive+Enrichment',
                    cases4print[5]:'Adaptive',
                    cases4print[6]:'A-ADM'}
        markers4print={cases4print[0]:'o',
                    cases4print[1]:'v',
                    cases4print[2]:'D',
                    cases4print[3]:'s',
                    cases4print[4]:'p',
                    cases4print[5]:'x',
                    cases4print[6]:'1'}

        markeverys=3
        markevery=300
        for case in all_cases:
            case_name=case[0]
            case_data=case[1]
            case_data['vpi']
            if case_name[3]!='i':
                style='-.'
            else:
                style='-'

            if variable!='vpi' and variable!='vpis_for_save' and (case_name in cases4print):
                if case_name=='finescale':
                    if variable!='n1_adm' and variable!='el2' and variable!='elinf' and variable!='es_L2' and variable!='es_Linf':
                        plt.plot(100*case_data['vpi'], case_data[variable],label=names4print[case_name], marker=markers4print[case_name],markevery=markevery, markersize=10)
                else:
                    if variable=='n1_adm':
                        plt.plot(100*case_data['vpi'], 100*case_data[variable][:-1]/case_data[variable][-1], style, label=names4print[case_name], marker=markers4print[case_name],markevery=markevery, markersize=10)
                        plt.xlabel(units['vpi'])
                        plt.ylabel(units[variable])
                    elif variable!='es_L2' and variable!='es_Linf':
                        if case_data[variable].min()<ymin:
                            ymin=case_data[variable].min()
                        if case_data[variable].max()>ymax:
                            ymax=case_data[variable].max()

                        plt.plot(100*case_data['vpi'], case_data[variable],style,label=names4print[case_name], marker=markers4print[case_name],markevery=markevery, markersize=10)
                        plt.xlabel(units['vpi'])
                        plt.ylabel(units[variable])
                    else: # Saturation variables
                        if case_data[variable].min()<ymin:
                            ymin=case_data[variable].min()
                        if case_data[variable].max()>ymax:
                            ymax=case_data[variable].max()
                        try:
                            plt.plot(100*case_data['vpis_for_save'], case_data[variable],style,label=names4print[case_name], marker=markers4print[case_name],markevery=markeverys, markersize=10)
                        except:
                            import pdb; pdb.set_trace()
                        plt.xlabel(units['vpis_for_save'])
                        plt.ylabel(units[variable])


        if variable!='vpi' and variable!='vpis_for_save':
            if (variable=='el2' or variable=='elinf' or variable=='es_L2' or variable=='es_Linf') and False:
                if ymin>0:
                    ymin=10**int(np.log10(ymin))
                if ymax>0:
                    ymax=10**int(np.log10(ymax)+1)
                # plt.yscale('log')
                y_ticks=plt.gca().get_yticks()[:-1]
                y_ticks=np.concatenate(np.vstack([y_ticks,np.append(y_ticks[1:]/2,y_ticks[-1])]).T)[:-1]
                plt.gca().set_yticks(y_ticks)
                ticks=plt.gca().get_yticks()

                pos_ymax=(ticks<ymax).sum()
                t_ymin=ticks[ticks>=ymin].min()
                # import pdb; pdb.set_trace()
                t_ymax=ticks[pos_ymax]

                plt.yticks(ticks)

                plt.gca().set_yticklabels(['{:.0f}%'.format(x) for x in plt.gca().get_yticks()])
                plt.ylim(t_ymin,t_ymax)

            plt.gca().tick_params(which='minor', length=10)
            plt.gca().tick_params(which='major', length=15)

            plt.grid()
            plt.legend()
            plt.gcf().set_size_inches(15,15)
            plt.savefig('results/'+folder+'/'+variable+'.svg')


def print_results_3(all_cases):
    units={'vpi':'PVI [%]','wor':'wor []','t_comp':'comp_time [s]','delta_t':'time-step []',
        'n1_adm':r'$N^{A-ADM}/N^f$[%]','el2':r'$||e_p||_2$'+' [%]','elinf':r'$||e_p||_\infty$'+' [%]', 'es_L2':r'$||e_s||_2$'+' [%]',
        'es_Linf':r'$||e_s||_\infty$'+' [%]', 'vpis_for_save':'PVI [%]','coupl':'Percentage of enrichment [%]',
        'refinement':'Percentage at fine-scale [%]', 'ep_haji_L2':'ep_rel_ms_L2[]',
        'ep_haji_Linf':'ep_rel_ms_Linf[]','er_L2':r'$||e_r||_2$'+' []','er_Linf':r'$||e_r||_\infty$'+' []',
        'ev_L2':r'$||e_v||_2$'+' [%]','ev_Linf':r'$||e_v||_\infty$'+' [%]'}

    variables=all_cases[1][1].keys()

    cases4print=['finescale',
                'neta_inf_alpha_1.0_type_uni_beta_3.0_delta_0.1_CR_0',
                'neta_inf_alpha_1.0_type_uni_beta_3.0_delta_0.3_CR_0',
                'neta_inf_alpha_1.0_type_uni_beta_3.0_delta_0.5_CR_0',
                'neta_inf_alpha_1.0_type_uni_beta_3.0_delta_1.1_CR_0',
                'neta_inf_alpha_1.0_type_uni_beta_3.0_delta_1.1_CR_0s',
                'neta_inf_alpha_1.0_type_uni_beta_inf_delta_1.1_CR_0s',
                'neta_inf_alpha_1.0_type_uni_beta_3.0_delta_1.1_CR_0s',
                'neta_inf_alpha_1.0_type_uni_beta_3.0_delta_1.1_CR_0s',
                'neta_inf_alpha_1.0_type_uni_beta_inf_delta_1.1_CR_0s',
                ]



    names4print={cases4print[0]:'Reference',
                cases4print[1]:'AMS   : '+r'$\gamma^{lim}=0.1$',
                cases4print[2]:'AMS   : '+r'$\gamma^{lim}=0.3$',
                cases4print[3]:'AMS   : '+r'$\gamma^{lim}=0.5$',
                cases4print[4]:'AMS   : '+r'$\gamma^{lim}=\infty$',
                cases4print[5]:'A-AMS: '+r'$\gamma^{lim}=0.1$',
                cases4print[6]:'A-AMS: '+r'$\gamma^{lim}=0.3$',
                cases4print[7]:'A-AMS: '+r'$\gamma^{lim}=0.5$',
                cases4print[8]:'A-AMS: '+r'$\gamma^{lim}=\infty$',
                }


    markers4print={cases4print[0]:'*',
                    cases4print[1]:'$0.1$',
                    cases4print[2]:'$0.3$',
                    cases4print[3]:'$0.5$',
                    cases4print[4]:'$\infty$',
                    cases4print[5]:'$0.1$',
                    cases4print[6]:'$0.3$',
                    cases4print[7]:'$0.5$',
                    cases4print[8]:'$\infty$',}


    colors4print={cases4print[0]:'magenta',
                    cases4print[1]:'green',
                    cases4print[2]:'black',
                    cases4print[3]:'red',
                    cases4print[4]:'blue',
                    cases4print[5]:'green',
                    cases4print[6]:'black',
                    cases4print[7]:'red',
                    cases4print[8]:'blue',}

    linestyles4print={cases4print[0]:'solid',
                    cases4print[1]:'solid',
                    cases4print[2]:'solid',
                    cases4print[3]:'solid',
                    cases4print[4]:'solid',
                    cases4print[5]:(0,(1,2)),
                    cases4print[6]:(0,(1,2,2,2)),
                    cases4print[7]:(0,(1,2,4,2)),
                    cases4print[8]:(0,(1,2,8,2)),
                    }

    markeverys=15
    markevery=1500
    markersize=40

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
        plt.rcParams['axes.edgecolor'] = 'black'
        plt.rcParams['axes.linewidth'] = 5
        plt.gca().tick_params(which='minor', length=10, width=5)
        plt.gca().tick_params(which='major', length=25, width=5)

        ymin=np.inf
        ymax=-np.inf

        max_neta=1
        for case in all_cases:
            case_name=case[0]
            case_data=case[1]
            if case_name!='finescale':
                try:
                    single_vars[variable].append(case_data[variable][0])
                except:
                    pass
                if variable=='el2':
                    input_params=case_name.split('_')
                    for i, param in enumerate(input_params):
                        if param in names_single_vars:

                            if input_params[i+1]=='inf':
                                single_vars[param].append(max_neta)
                            else:
                                single_vars[param].append(float(input_params[i+1]))
                style=''
            else:
                style='-.'
            if variable!='vpi' and variable!='vpis_for_save' and (case_name in cases4print):

                if case_name=='finescale':
                    if variable!='n1_adm' and variable!='el2' and variable!='elinf' and variable!='es_L2'and variable!='es_Linf' and variable!='refinement' and variable!='coupl':

                        try:
                            plt.plot(100*case_data['vpi'], case_data[variable],label=names4print[case_name], linewidth=5, marker=markers4print[case_name],
                             markevery=markevery,color=colors4print[case_name], markersize=markersize, linestyle=linestyles4print[case_name])
                        except:
                            pass
                else:
                    if variable=='n1_adm':
                        plt.plot(100*case_data['vpi'], 100*case_data[variable][:-1]/case_data[variable][-1], style, label=names4print[case_name], linewidth=5,
                         marker=markers4print[case_name], markevery=markevery, color=colors4print[case_name], markersize=markersize, linestyle=linestyles4print[case_name])
                        plt.xlabel(units['vpi'])
                        plt.ylabel(units[variable])

                    elif variable!='es_L2' and variable!='es_Linf':

                        if case_data[variable].min()<ymin:
                            ymin=case_data[variable].min()
                        if case_data[variable].max()>ymax:
                            ymax=case_data[variable].max()
                        try:
                            plt.plot(100*case_data['vpi'], case_data[variable],style,label=names4print[case_name], linewidth=5, marker=markers4print[case_name],
                             markevery=markevery, color=colors4print[case_name], markersize=markersize, linestyle=linestyles4print[case_name])
                        except:
                            pass

                        plt.xlabel(units['vpi'])
                        plt.ylabel(units[variable])

                    else: # Saturation variables
                        if case_data[variable].min()<ymin:
                            ymin=case_data[variable].min()
                        if case_data[variable].max()>ymax:
                            ymax=case_data[variable].max()

                        plt.plot(100*case_data['vpis_for_save'], case_data[variable],style,label=names4print[case_name], linewidth=5,markevery=markeverys,
                        marker=markers4print[case_name], color=colors4print[case_name], markersize=markersize, linestyle=linestyles4print[case_name])

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

                y_ticks=plt.gca().get_yticks()[:-1]
                y_ticks=np.concatenate(np.vstack([y_ticks,np.append(y_ticks[1:]/2,y_ticks[-1])]).T)[:-1]
                plt.gca().set_yticks(y_ticks)
                ticks=plt.gca().get_yticks()

                pos_ymax=(ticks<ymax).sum()

            plt.grid(linewidth=3)

            plt.legend()
            plt.gcf().set_size_inches(20,20)
            plt.tick_params(colors='black')
            ss=plt.gca().get_legend_handles_labels()
            '''Organizes plot legend'''
            handles, labels = plt.gca().get_legend_handles_labels()
            ind_sort=np.argsort(labels)
            slabels=np.array(labels)[ind_sort].tolist()
            shandles=np.array(handles)[ind_sort].tolist()
            plt.legend(shandles, slabels)
            '''end'''
            plt.savefig('results/'+folder+'/'+variable+'.svg', bbox_inches='tight')

    for var in single_vars:
        abcissa_var='alpha'
        control_variable='beta'

        all_abcissa=np.array(single_vars[abcissa_var])
        pos=all_abcissa<10000
        all_abcissa=all_abcissa[pos]
        all_ordenada=np.array(single_vars[var])[pos]
        control_parameters=single_vars[control_variable]

        ordenadas=[]
        abcissas=[]
        control_vals=np.unique(control_parameters)

        for p in np.unique(control_parameters):
            orden=all_ordenada[control_parameters==p]
            abcis=all_abcissa[control_parameters==p]
            if abcis.max()<10000:
                if (all_abcissa==max_neta).sum()>0:
                    orden=np.concatenate([[all_ordenada[all_abcissa==max_neta][0]],orden])

                abcis=np.concatenate([[max_neta],abcis])
            ordenadas.append(orden)
            abcissas.append(abcis)

        plt.close('all')

        linear_yaxis=['refinement','coupl', 'ev_L2','ev_Linf']
        if var in ['elinf','el2', 'ev_L2', 'ev_Linf', 'refinement', 'coupl']:

            plt.close('all')
            fig=plt.figure()
            plt.grid()
            plt.xlabel(r'$\epsilon$ []',labelpad=-20)
            # if var not in ['neta', 'alpha', 'beta', 'delta']:
            #     plt.ylabel(units[var])#, labelpad=-20)

            control=0

            for abcissa, ordenada in zip(abcissas, ordenadas):
                plt.gca().tick_params(which='minor', length=10)
                plt.gca().tick_params(which='major', length=15)
                ind_sort=np.argsort(abcissa)
                abcissa=abcissa[ind_sort]
                try:
                    ordenada=ordenada[ind_sort]
                except:
                    ordenada=np.repeat(ordenada[0],len(ind_sort))

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

                plt.plot(abc, ord,label=label.format(control_val),marker='o',markersize=15, linewidth=5)

                plt.rcParams['axes.edgecolor'] = 'black'
                plt.rcParams['axes.linewidth'] = 1

                positions = np.unique(all_abcissa.astype(int))
                positions=positions[positions!=30]
                labels = positions.astype(str)


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

                if var in linear_yaxis:
                    pp.append(int(all_ordenada.min()))
                    pp.append(int(all_ordenada.max()))


                close_values=[]
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

                control+=1
            # plt.xscale('log')
            plt.gca().tick_params(which='minor', length=10)
            plt.gca().tick_params(which='major', length=15)

            plt.grid(lw=3)
            # plt.legend()

            plt.gcf().set_size_inches(15,15)

            plt.savefig('results/single_phase/'+var+'.svg', bbox_inches='tight')

def collect_single_phase_data(all_cases):
    units={'vpi':'PVI [%]','wor':'wor []','t_comp':'comp_time [s]','delta_t':'time-step []',
        'n1_adm':r'$N^{A-ADM}/N^f$[%]','el2':r'$||e_p||_2$'+' [%]','elinf':r'$||e_p||_\infty$'+' [%]', 'es_L2':r'$||e_s||_2$'+' [%]',
        'es_Linf':r'$||e_s||_\infty$'+' [%]', 'vpis_for_save':'PVI [%]','coupl':'Percentage of enrichment [%]',
        'refinement':'Percentage at fine-scale [%]', 'ep_haji_L2':'ep_rel_ms_L2[]',
        'ep_haji_Linf':'ep_rel_ms_Linf[]','er_L2':r'$||e_r||_2$'+' []','er_Linf':r'$||e_r||_\infty$'+' []',
        'ev_L2':r'$||e_v||_2$'+' [%]','ev_Linf':r'$||e_v||_\infty$'+' [%]'}
    variables=all_cases[1][1].keys()
    single_vars={}
    for u in units.keys():
        single_vars[u]=[]
    single_vars['neta']=[]
    single_vars['alpha']=[]
    single_vars['beta']=[]
    single_vars['delta']=[]
    single_vars['CR']=[]
    names_single_vars=['neta', 'beta', 'alpha', 'delta', 'CR']

    for variable in single_vars.keys():
        for case in all_cases:
            case_name=case[0]
            case_data=case[1]
            if case_name!='finescale' and variable in names_single_vars:
                input_params=np.array(case_name.split('_'))
                pos=np.where(input_params==variable)[0][0]+1
                single_vars[variable].append(float(input_params[pos]))
            elif case_name!='finescale':
                try:
                    single_vars[variable].append(case_data[variable][0])
                except:
                    pass
        single_vars[variable]=np.array(single_vars[variable])
    return single_vars

def round_to_1(x):

    vals=-np.floor(np.log10(abs(x))).astype('int')
    for i in range(len(x)):
        x[i]=np.round(x[i], vals[i])
    return x

def format_plot(scales, abcissas, ordenadas):
    x_scale, y_scale = scales.split('_')
    if x_scale=='lin':
        x_scale='linear'
    if y_scale=='lin':
        y_scale='linear'
    plt.xscale(x_scale)
    plt.yscale(y_scale)

    plt.grid(which='major', lw=2, color='black')
    plt.grid(which='minor', lw=1, color='gray')

    if y_scale=='log':
        major_ticks=np.log10(ordenadas).astype('int')
        major_ticks=10**np.arange(major_ticks.min(),major_ticks.max()+1)

        plt.gca().yaxis.set_major_locator(ticker.FixedLocator(major_ticks))
        plt.gca().set_yticklabels(['{:.0f}%'.format(x) for x in np.concatenate([major_ticks])])

        minor_ticks=np.unique(np.array(ordenadas).astype('int'))
        s_minor_ticks=np.arange(major_ticks.max(),max(2.01*major_ticks.max(),ordenadas.max()*1.101),10**(np.log10(major_ticks.max())))[1:]
        minor_ticks=round_to_1(minor_ticks)
        a_minor_ticks=np.concatenate([np.arange(major_ticks[i],major_ticks[i+1],major_ticks[i])[1:] for i in range(len(major_ticks)-1)])
        a_minor_ticks=np.concatenate([a_minor_ticks, s_minor_ticks])

        amt=a_minor_ticks.astype('str')
        minor_ticks=np.append(minor_ticks,a_minor_ticks.max())
        a_minor_formats=[]
        plot_mantissa=['2', '4', '7']
        plt.yticks(fontsize=50)
        for i in range(len(a_minor_ticks)):
            if (amt[i][0] in plot_mantissa) or (a_minor_ticks[i]==a_minor_ticks.max()):
                a_minor_formats.append('{:.0f}%'.format(a_minor_ticks[i]))
            else:
                a_minor_formats.append('')

                plt.gca().yaxis.set_minor_locator(ticker.FixedLocator(a_minor_ticks))
                plt.gca().yaxis.set_minor_formatter(ticker.FixedFormatter(a_minor_formats))
                plt.ylim(minor_ticks.min(),minor_ticks.max()*1.1)

    plt.gcf().set_size_inches(15,15)


def print_results_single(all_cases):
    single_vars = collect_single_phase_data(all_cases)

    plot_vars=[[    'neta',    'neta'],     # Abcissas
               [     'el2',   'elinf'],     # Ordenadas
               [ 'lin_lin', 'lin_log']]     # Escalas dos Eixos

    for abcissa_name, ordenada_name, scales in zip(plot_vars[0],plot_vars[1], plot_vars[2]):
        plt.close('all')
        abcissas=single_vars[abcissa_name]
        ordenadas=single_vars[ordenada_name]

        ind_sort=np.argsort(abcissas)
        abcissas=abcissas[ind_sort]
        ordenadas=ordenadas[ind_sort]
        # import pdb; pdb.set_trace()
        plt.plot(abcissas, ordenadas, lw=5)
        format_plot(scales, abcissas, ordenadas)
        plt.savefig('results/single_phase/'+ordenada_name+'.svg', bbox_inches='tight')
