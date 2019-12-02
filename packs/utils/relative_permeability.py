from .. import directories as direc
import numpy as np

def brooks_and_corey(saturations):

    def _stemp(S, Swc, Sor):
        return (S - Swc) / (1 - Swc - Sor)

    def _krw(S_temp, krw0, n_w):
        return krw0*(np.power(S_temp, n_w))

    def _kro(S_temp, kro0, n_o):
        return kro0*(np.power(1 - S_temp, n_o))

    Sor = float(direc.data_loaded['biphasic_data']['Sor'])
    Swc = float(direc.data_loaded['biphasic_data']['Swc'])
    n_w = float(direc.data_loaded['biphasic_data']['n_w'])
    n_o = float(direc.data_loaded['biphasic_data']['n_o'])
    krw0 = float(direc.data_loaded['biphasic_data']['krw0'])
    kro0 = float(direc.data_loaded['biphasic_data']['kro0'])
    n = len(saturations)
    ids = np.arange(n)
    ids_fora = ids[(saturations < 0) | (saturations > 1)]
    if len(ids_fora) > 0:
        raise ValueError('saturacao errada')

    krw = np.zeros(n)
    kro = krw.copy()
    ids1 = ids[saturations > (1 - Sor)]
    ids2 = ids[saturations < Swc]
    ids_r = np.setdiff1d(ids, np.union1d(ids1, ids2))

    krw[ids1] = np.repeat(krw0, len(ids1))
    kro[ids2] = np.repeat(kro0, len(ids2))
    Stemp = _stemp(saturations[ids_r], Swc, Sor)
    krw[ids_r] = _krw(Stemp, krw0, n_w)
    kro[ids_r] = _kro(Stemp, kro0, n_o)

    return krw, kro
