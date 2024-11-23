import numpy as np


class GlobalDebug:
    _global_debug = False
    _debug = False

    @classmethod
    def return_ones(cls, vec: np.ndarray):
        k = 1
        if cls._debug is True or cls._global_debug is True:
            return k * np.ones(vec.shape)
        else:
            return vec
