import numpy as np
import os
import matplotlib.pyplot as plt
import shutil
import matplotlib
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
    os.system("python testting2_biphasic.py")
    neta_lim_dual_values=     [ np.inf,    1.0,    5.0,   10.0,   50.0,  100.0]#,    2.0,    5.0,   10.0,  100.0,   500.0,1000.0, np.inf]
    neta_lim_finescale_values=[ np.inf, np.inf, np.inf, np.inf, np.inf, np.inf]#, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf]
    type_of_refinement_values=[  'uni',  'uni',  'uni',  'uni',  'uni',  'uni']#,  'uni',  'uni',  'uni',  'uni',  'uni',  'uni',  'uni']
    phiK_raz_lim_values=      [ np.inf, np.inf, np.inf, np.inf, np.inf, np.inf]#, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf]
    delta_sat_max=            [    1.1,    1.1,    1.1,    1.1,    1.1,    1.1]#,    1.1,    1.1,    1.1,    1.1,    1.0,    1.0,    1.0]
    for i in range(len(neta_lim_dual_values)):
        np.save('flying/delta_sat_max.npy',np.array([delta_sat_max[i]]))
        np.save('flying/neta_lim_finescale.npy',np.array([neta_lim_finescale_values[i]]))
        np.save('flying/neta_lim_dual.npy',np.array([neta_lim_dual_values[i]]))
        np.save('flying/type_of_refinement.npy',np.array([type_of_refinement_values[i]]))
        np.save('flying/phiK_raz_lim.npy',np.array([phiK_raz_lim_values[i]]))
        ms_case='neta_'+str(neta_lim_dual_values[i])+'_alpha_'+str(neta_lim_finescale_values[i])+'_'+type_of_refinement_values[i]+'_beta_'+str(phiK_raz_lim_values[i])+'_delta_'+str(delta_sat_max[i])+'/'
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

def get_axis_lims(ordenada,max_raz):
    if (ordenada==0).sum()==0:
        reasons=ordenada[:-1]/ordenada[1:]
        reasons[1/reasons>reasons]=1/reasons[1/reasons>reasons]
        reasons=np.concatenate([[1],reasons])
        ind_superates=np.arange(len(ordenada))[(reasons>max_raz) | (reasons<1/max_raz)]
        i_prev=0
        ylims=[]

        for ind in ind_superates:
            vals=np.array([ordenada[i_prev],ordenada[ind-1]])
            if vals.max()<100:
                ylims.append((vals.min()*0.99,vals.max()*1.01))
            else:
                ylims.append((vals.min()-0.99,vals.max()+1.01))
            i_prev=ind
        # if i_prev<=len(ordenada)-1:
        if ordenada[-1]<100:
            ylims.append((ordenada[-1]*0.99,ordenada[-1]*1.01))
        else:
            ylims.append((ordenada[-1]-0.99,ordenada[-1]+1.01))
    else:
        if ordenada.min()!=ordenada.max():
            ylims=[(ordenada.min(),ordenada.max())]
        else:
            ylims=[(0,1)]
    # print(ylims)
    return(ylims)

def print_results(all_cases):
    units={'vpi':'vpi [%]','wor':'wor []','t_comp':'comp_time [s]','delta_t':'time-step []',
        'n1_adm':'Nadm/Nf[%]','el2':r'$||e_p||_2$'+' [%]','elinf':r'$||e_p||_\infty$'+' [%]', 'es_L2':'es_L2 [%]',
        'es_Linf':'es_Linf [%]', 'vpis_for_save':'vpi [%]','coupl':"Percentage of enrichment [%]",
        'refinement':"Percentage at fine-scale [%]", 'ep_haji_L2':'ep_rel_ms_L2[]',
        'ep_haji_Linf':'ep_rel_ms_Linf[]','er_L2':r'$||e_r||_2$'+' []','er_Linf':r'$||e_r||_\infty$'+' []'}
    variables=all_cases[0][1].keys()
    single_vars={}
    for u in units.keys():
        single_vars[u]=[]
    single_vars['neta']=[]
    for variable in variables:
        plt.close('all')
        ymin=np.inf
        ymax=-np.inf

        for case in all_cases:
            case_name=case[0]
            case_data=case[1]


            if case_name!='finescale':
                single_vars[variable].append(case_data[variable][0])
                if case_name[5:8]!='inf':
                    if variable=='el2':
                        single_vars['neta'].append(float(case_name[5:8]))
                else:
                    if variable=='el2':
                        single_vars['neta'].append(1000.0)

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
                        try:
                            plt.plot(100*case_data['vpi'], case_data[variable],style,label=case_name)
                        except:
                            import pdb; pdb.set_trace()
                        plt.xlabel(units['vpi'])
                        plt.ylabel(units[variable])

                    else: # Saturation variables
                        if case_data[variable].min()<ymin:
                            ymin=case_data[variable].min()
                        if case_data[variable].max()>ymax:
                            ymax=case_data[variable].max()
                        # try:
                        plt.plot(100*case_data['vpis_for_save'], case_data[variable],style,label=case_name)
                        # except:
                            # import pdb; pdb.set_trace()
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
            plt.gcf().set_size_inches(20,10)
            plt.savefig('results/biphasic/'+variable+'.png')

    for var in single_vars:
        if var!='neta':
            plt.close('all')
            abcissa=np.array(single_vars['neta'])
            pos=abcissa<10000
            abcissa=abcissa[pos]
            ordenada=np.array(single_vars[var])[pos]
            ind_sort=np.argsort(abcissa)
            abcissa=abcissa[ind_sort]
            ordenada=ordenada[ind_sort]
            if 'haji' in var:
                ordenada[abcissa==1000]=1

            xlims=get_axis_lims(abcissa,9)
            ylims=get_axis_lims(ordenada,10)
            ylims[0]=(0,ylims[0][1])
            xlims[0]=(0,xlims[0][1])

            fig=plt.figure()
            bax=brokenaxes(xlims=xlims, ylims=ylims)

            # import pdb; pdb.set_trace()
            # if ordenada[-2]!=0 and ordenada[-1]/ordenada[-2]>10 and var!='coupl':
            #
            #     ylims=((0,ordenada[-2]*1.01),(ordenada[-1]*0.99,ordenada[-1]*1.01))
            #
            #     bax=brokenaxes(xlims=[(0,abcissa[-2]*1.01),(999,1001)], ylims=ylims)
            #
            # else:
            #     bax=brokenaxes(xlims=((0,abcissa[-2]*1.01),(999,1001)))
            # if len(ylims1)>1:
            #     import pdb; pdb.set_trace()
            # bax=brokenaxes(xlims=[(0,abcissa[-2]*1.01),(999,1001)], ylims=ylims)

            plt.gca().tick_params(which='minor', length=10)
            plt.gca().tick_params(which='major', length=15)
            plt.gcf().set_size_inches(20,15)
            bax.set_xlabel(r'$\epsilon$ []',labelpad=40)
            bax.set_ylabel(units[var], labelpad=70)
            bax.plot(abcissa, ordenada)
            bax.scatter(abcissa, ordenada, s=60)
            bax.grid()
            plt.savefig('results/single_phase/'+var+'.png')

    import pdb; pdb.set_trace()
