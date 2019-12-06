from ..flux_schemes.tpfa_scheme import TpfaScheme

class FineScaleTpfaPressureSolver(TpfaScheme):

    def __init__(self, data_impress, elements_lv0, data_name: str='FineScaleTpfaPressureSolver.npz', load=False):
        super().__init__(data_impress, elements_lv0, data_name=data_name, load=load)
