import numpy as np
import os
import matplotlib.pyplot as plt
import shutil
import matplotlib
font = {'size'   : 22}
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
    vpis_for_save=np.arange(0.0,0.501,0.01)

    np.save('flying/vpis_for_save.npy',vpis_for_save)
    os.system("python testting2_biphasic.py")
    neta_lim_dual_values=     [ np.inf, np.inf, np.inf, np.inf,    1.0,    1.0,    1.0,    1.0,    1.0, np.inf, np.inf]
    neta_lim_finescale_values=[ np.inf, np.inf,    1.0,    1.0, np.inf, np.inf,    1.0,    1.0,    1.0,    1.0,    1.0]
    type_of_refinement_values=['n_uni',  'uni',  'uni',  'uni',  'uni',  'uni',  'uni',  'uni',  'uni',  'uni',  'uni']
    phiK_raz_lim_values=      [ np.inf, np.inf, np.inf,      3, np.inf,      3, np.inf,      3,      3, np.inf,    3.0]
    delta_sat_max=            [    0.1,    0.1,    0.1,    0.1,    0.1,    0.1,    0.1,    0.1,    2.0,    2.0,    2.0]
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

def print_results(all_cases):

    units={'vpi':'vpi[%]','wor':'wor[%]','t_comp':'comp_time[s]','delta_t':'time-step[]', 'n1_adm':'Nadm/Nf[%]','el2':'ep_L2[%]','elinf':'ep_Linf[%]', 'es_L2':'es_L2[%]', 'es_Linf':'es_Linf[%]', 'vpis_for_save':'vpi[%]'}
    variables=all_cases[0][1].keys()

    for variable in variables:
        plt.close('all')
        ymin=np.inf
        ymax=-np.inf

        for case in all_cases:
            case_name=case[0]
            case_data=case[1]
            # case_data['vpi']
            if case_name[5]!='i':
                style='-.'
            else:
                style='-'
            if variable!="vpi" and variable!='vpis_for_save':
                if case_name=='finescale':
                    if variable!="n1_adm" and variable!="el2" and variable!="elinf" and variable!='es_L2' and variable!='es_Linf':
                        plt.plot(100*case_data['vpi'], case_data[variable],label=case_name)
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
                        try:
                            plt.plot(100*case_data['vpis_for_save'], case_data[variable],style,label=case_name)
                        except:
                            import pdb; pdb.set_trace()
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
                # import pdb; pdb.set_trace()
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
