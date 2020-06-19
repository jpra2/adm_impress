import numpy as np
import os
import matplotlib.pyplot as plt
import shutil
import matplotlib

def remove_previous_files():
    # Cleans the previous results
    shutil.rmtree('results/biphasic/ms')
    shutil.rmtree('results/biphasic/finescale')
    # creates new folders
    os.mkdir('results/biphasic/ms')
    os.mkdir('results/biphasic/finescale')
    os.mkdir('results/biphasic/finescale/vtks')

def run_test_cases():
    # neta_lim_dual_values=[np.inf, np.inf, np.inf,1.0, 1.0, 1.0, np.inf, 1.0]
    # neta_lim_finescale_values=[np.inf, 5.0, 0.5, np.inf, 5.0, 0.5, np.inf, np.inf]
    # type_of_refinement_values=['uni','uni','uni', 'uni', 'uni', 'uni', 'n_uni', 'n_uni']
    # phiK_raz_lim_values=[3, 3, 3, 3, 3, 3, 3, 3]
    neta_lim_dual_values=[np.inf, np.inf,np.inf, 1.0, 1.0,1.0]
    neta_lim_finescale_values=[20.0, 100.0,1000.0, 20.0, 100.0, 1000.0]
    type_of_refinement_values=['uni','uni','uni', 'uni','uni', 'uni', 'uni']
    phiK_raz_lim_values=[5, 5, 5, 5, 5, 5]
    for i in range(len(neta_lim_dual_values)):
        np.save('flying/neta_lim_finescale.npy',np.array([neta_lim_finescale_values[i]]))
        np.save('flying/neta_lim_dual.npy',np.array([neta_lim_dual_values[i]]))
        np.save('flying/type_of_refinement.npy',np.array([type_of_refinement_values[i]]))
        np.save('flying/phiK_raz_lim.npy',np.array([phiK_raz_lim_values[i]]))
        ms_case='n_d'+str(neta_lim_dual_values[i])+'_nf_'+str(neta_lim_finescale_values[i])+'_'+type_of_refinement_values[i]+'_sens_'+str(phiK_raz_lim_values[i])+'/'
        os.makedirs('results/biphasic/ms/'+ms_case,exist_ok=True)
        os.makedirs('results/biphasic/ms/'+ms_case+'vtks',exist_ok=True)
        np.save('flying/ms_case.npy',np.array([ms_case]))
        os.system("python ADM_b.py")

    os.system("python testting2_biphasic.py")

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
    font = {'size'   : 22}
    units={'vpi':'vpi[%]','wor':'wor[%]','t_comp':'comp_time[s]','delta_t':'time-step[]', 'n1_adm':'Nadm/Nf[%]','el2':'ep_L2[%]','elinf':'ep_Linf[%]'}
    variables=all_cases[0][1].keys()
    ymin=np.inf
    ymax=-np.inf
    for variable in variables:
        plt.close('all')
        for case in all_cases:
            case_name=case[0]
            case_data=case[1]
            case_data['vpi']
            if case_name[3]!='i':
                style='-.'
            else:
                style='-'
            if variable!="vpi":
                if case_name=='finescale':
                    if variable!="n1_adm" and variable!="el2" and variable!="elinf":
                        plt.plot(100*case_data['vpi'], case_data[variable],label=case_name)
                else:
                    if variable=='n1_adm':
                        plt.plot(100*case_data['vpi'], 100*case_data[variable][:-1]/case_data[variable][-1], style, label=case_name)
                        plt.xlabel(units['vpi'])
                        plt.ylabel(units[variable])
                    else:

                        if case_data[variable].min()<ymin:
                            ymin=case_data[variable].min()
                        elif case_data[variable].max()>ymax:
                            ymax=case_data[variable].max()

                        plt.plot(100*case_data['vpi'], case_data[variable],style,label=case_name)
                        plt.xlabel(units['vpi'])
                        plt.ylabel(units[variable])


        if variable!='vpi':
            if variable=='el2' or variable=='elinf':
                plt.yscale('log')
                y_ticks=plt.gca().get_yticks()[:-1]
                y_ticks=np.concatenate(np.vstack([y_ticks,np.append(y_ticks[1:]/2,y_ticks[-1])]).T)[:-1]
                plt.gca().set_yticks(y_ticks)
                ticks=plt.gca().get_yticks()
                # import pdb; pdb.set_trace()
                t_ymin=ticks[ticks>=ymin].min()*10
                t_ymax=ticks[ticks<=ymax].max()*10

                # import pdb; pdb.set_trace()
                plt.yticks(ticks)
                plt.gca().tick_params(which='minor', length=4, color='r')                
                plt.gca().set_yticklabels(['{:.0f}%'.format(x) for x in plt.gca().get_yticks()])
                plt.ylim((t_ymin,t_ymax))



            matplotlib.rc('font', **font)
            plt.grid()
            plt.legend()
            plt.gcf().set_size_inches(20,10)
            plt.savefig('results/biphasic/'+variable+'.png')
