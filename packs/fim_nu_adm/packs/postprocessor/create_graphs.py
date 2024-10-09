import matplotlib.pyplot as plt
import numpy as np
case='label_85'

def print_results():
    p_ms=np.load('results/'+case+'/pressures_multilevel.npy')
    p_fs=np.load('results/'+case+'/pressures_finescale.npy')
    s_ms=np.load('results/'+case+'/saturations_multilevel.npy')
    s_fs=np.load('results/'+case+'/saturations_finescale.npy')
    n1_adm=np.load('results/'+case+'/n1_adm_multilevel.npy')
    phis=np.load('results/'+case+'/porosities.npy')
    p_ms_cpr=np.load('results/'+case+'/pressures_multilevel_CPR.npy')
    n1_adm_cpr=np.load('results/'+case+'/n1_adm_multilevel_CPR.npy')
    s_ms_cpr=np.load('results/'+case+'/saturations_multilevel_CPR.npy')
    n1_adm_cpr=np.load('results/'+case+'/n1_adm_multilevel_CPR.npy')
    ep_l2_CPR=np.load('results/'+case+'/ep_l2_multilevel_CPR.npy')
    ep_li_CPR=np.load('results/'+case+'/ep_li_multilevel_CPR.npy')
    es_l1_CPR=np.load('results/'+case+'/es_l1_multilevel_CPR.npy')
    # ep_l2=np.load('results/'+case+'/ep_l2_multilevel.npy')
    # ep_li=np.load('results/'+case+'/ep_li_multilevel.npy')
    # es_l1=np.load('results/'+case+'/es_l1_multilevel.npy')
    # import pdb; pdb.set_trace()
    # import pdb; pdb.set_trace()
    nf=len(p_ms[0])
    vpi_fs=100*(s_fs*phis).sum(axis=1)/phis.sum()
    vpi_ms=100*(s_ms*phis).sum(axis=1)/phis.sum()
    vpi_ms_cpr=100*(s_ms_cpr*phis).sum(axis=1)/phis.sum()
    pva=100*n1_adm/nf
    pva_cpr=100*n1_adm_cpr/nf
    # import pdb; pdb.set_trace()
    if len(vpi_fs)>len(vpi_ms):
        inds=np.searchsorted(vpi_fs,vpi_ms)-1
        vpi_fs=vpi_fs[inds]
        p_fs=p_fs[inds]
        s_fs=s_fs[inds]
    else:
        inds=np.searchsorted(vpi_ms,vpi_fs)-1
        vpi_ms=vpi_ms[inds]
        p_ms=p_ms[inds]
        s_ms=s_ms[inds]
        pva=pva[inds]

    if len(vpi_fs)>len(vpi_ms_cpr):
        inds=np.searchsorted(vpi_fs,vpi_ms_cpr)-1
        vpi_fs_cpr=vpi_fs[inds]
        p_fs=p_fs[inds]
        s_fs=s_fs[inds]
        ep_l2_CPR=ep_l2_CPR[inds]
        ep_li_CPR=ep_li_CPR[inds]
        es_l1_CPR=es_l1_CPR[inds]
    else:
        inds=np.searchsorted(vpi_ms_cpr,vpi_fs)-1
        vpi_ms_cpr=vpi_ms_cpr[inds]
        p_ms_cpr=p_ms_cpr[inds]
        s_ms_cpr=s_ms_cpr[inds]
        pva_cpr=pva_cpr[inds]
        ep_l2_CPR=ep_l2_CPR[inds]
        ep_li_CPR=ep_li_CPR[inds]
        es_l1_CPR=es_l1_CPR[inds]
    # import pdb; pdb.set_trace()
    ep2=100*np.linalg.norm(p_fs-p_ms,axis=1)/np.linalg.norm(p_fs,axis=1)
    ep_inf=100*abs(p_fs-p_ms[0:len(p_fs)]).max(axis=1)/p_fs.max(axis=1)
    es_1=100*abs(s_fs-s_ms).sum(axis=1)/nf
    ep2_cpr=100*np.linalg.norm(p_fs-p_ms_cpr,axis=1)/np.linalg.norm(p_fs,axis=1)
    ep_inf_cpr=100*abs(p_fs-p_ms_cpr[0:len(p_fs)]).max(axis=1)/p_fs.max(axis=1)
    es_1_cpr=100*abs(s_fs-s_ms_cpr).sum(axis=1)/nf

    viz=np.concatenate(np.load('results/'+case+'/viz.npy'))
    # import pdb; pdb.set_trace()
    sp_fs=(s_fs[:,viz[0]]+s_fs[:,viz[1]])/2
    wor_fs=sp_fs/(1-sp_fs)
    sp_ms=(s_ms[:,viz[0]]+s_ms[:,viz[1]])/2
    wor_ms=sp_ms/(1-sp_ms)

    sp_ms_cpr=(s_ms_cpr[:,viz[0]]+s_ms_cpr[:,viz[1]])/2
    wor_ms_cpr=sp_ms_cpr/(1-sp_ms_cpr)
    # import pdb; pdb.set_trace()
    save_image([vpi_fs,vpi_ms,vpi_ms_cpr],[wor_fs,wor_ms,wor_ms_cpr],['Reference','FIM_NU_ADM', 'CPR_FIM_NU_ADM'],'wor []')#, marker=['.','.','*'])
    save_image([vpi_ms,vpi_ms_cpr],[ep2,ep2_cpr],['FIM_NU_ADM',"CPR_FIM_NU_ADM"],'ep2')
    save_image([vpi_ms,vpi_ms_cpr],[ep_inf,ep_inf_cpr],['FIM_NU_ADM',"CPR_FIM_NU_ADM"],'epinf')
    save_image([vpi_ms,vpi_ms_cpr],[es_1,es_1_cpr],['FIM_NU_ADM',"CPR_FIM_NU_ADM"],'es1')
    save_image([vpi_ms,vpi_ms_cpr],[pva,pva_cpr],['FIM_NU_ADM',"CPR_FIM_NU_ADM"],'pva')
    # import pdb; pdb.set_trace()
    # import pdb; pdb.set_trace()
    save_image([vpi_ms_cpr],[ep_l2_CPR],["CPR_FIM_NU_ADM"],'e2p')
    save_image([vpi_ms_cpr],[ep_li_CPR],["CPR_FIM_NU_ADM"],'eip')
    save_image([vpi_ms_cpr],[es_l1_CPR],["CPR_FIM_NU_ADM"],'e1s')
    import pdb; pdb.set_trace()
    # save_image([vpi_ms,vpi_ms_cpr],[ep_l2,ep_l2_CPR],['FIM_NU_ADM',"CPR_FIM_NU_ADM"],'ep_l2')
    # save_image([vpi_ms,vpi_ms_cpr],[ep_li,ep_li_CPR],['FIM_NU_ADM',"CPR_FIM_NU_ADM"],'ep_li')
    # save_image([vpi_ms,vpi_ms_cpr],[es_l1,es_l1_CPR],['FIM_NU_ADM',"CPR_FIM_NU_ADM"],'es_l1')



def save_image(abs,ords,labels,name):
    plt.close('all')
    marker="+"
    for ab,ord,label in zip(abs,ords,labels):
        if name[0:3]=='ep2' or name[0:3]=='e2p':
            ylabel=r'$||e_p||_2$ [%]'
        elif name[0:3]=='epi' or name[0:3]=='eip':
            ylabel=r'$||e_p||_\infty$ [%]'
        elif name[0:3]=='es1' or name[0:3]=='e1s':
            ylabel=r'$||e_s||_1$ [%]'
        elif name[0:3]=='pva':
            ylabel=r'$N^{NU-ADM}/N^f$ [%]'
        else:
            ylabel=name
        if name[0:3]=='wor':
            if label== 'CPR_FIM_NU_ADM':
                marker='*'
        if name[0]=='e':
            yscale='log'
        else:
            yscale='linear'
        # import pdb; pdb.set_trace()
        plt.plot(ab[1:],ord[1:],label=label,marker=marker)
        plt.yscale(yscale)
        # import pdb; pdb.set_trace()
        #[np.format_float_scientific(o,precision=3) for o in ord[1:]]
        marker='+'
        # plt.yscale('log')
    plt.legend()
    plt.xlabel('PVI [%]')

    plt.ylabel(ylabel)
    plt.grid()
    plt.savefig('results/'+case+'/'+name[0:3]+'.svg',transparent=True)
    plt.close('all')
