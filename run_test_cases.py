import numpy as np
import os
import matplotlib.pyplot as plt

neta_lim_dual_values=[np.inf, np.inf]
neta_lim_finescale_values=[0.05, 0.1]
for i in range(len(neta_lim_dual_values)):
    np.save('flying/neta_lim_finescale.npy',np.array([neta_lim_finescale_values[i]]))
    np.save('flying/neta_lim_dual.npy',np.array([neta_lim_dual_values[i]]))
    ms_case='n_d'+str(neta_lim_dual_values[i])+'_nf_'+str(neta_lim_finescale_values[i])+'/'
    os.makedirs('results/biphasic/ms/'+ms_case,exist_ok=True)
    np.save('flying/ms_case.npy',np.array([ms_case]))
#     os.system("python ADM_b.py")
#
# os.system("python testting2_biphasic.py")
cases_ms=[]
for root, dirs, files in os.walk("results/biphasic/ms"):
    for dir in dirs:
        variables={}
        for r, ds, fs in os.walk(os.path.join(root,dir)): #print(fs)
            for file in fs:
                var_name=os.path.splitext(file)[0]
                variables[var_name]=np.load(os.path.join(r,file))
        cases_ms.append([dir,variables])

finescale={}
for r, ds, fs in os.walk("results/biphasic/finescale"): #print(fs)
    for file in fs:
        var_name=os.path.splitext(file)[0]
        finescale[var_name]=np.load(os.path.join(r,file))

case_finescale=[['finescale',finescale]]
all_cases=np.array(cases_ms+case_finescale)
variables=all_cases[0][1].keys()
import pdb; pdb.set_trace()
for variable in variables:
    plt.close('all')
    for case in all_cases:
        case_name=case[0]
        case_data=case[1]
        if ((case_name!='finescale') and (variable!="n1_adm")) or (variable!="vpi"):
            plt.plot(case_data['vpi'], case_data[variable])

    plt.savefig('results/biphasic/'+variable+'.png')
