from packs.data_class.accumulative_array_data_manager import AccumulativeArrayDataManager

class CumulativeCompositionalDataArray(AccumulativeArrayDataManager):
    _dty = [('dvpi', float), ('vpi', float), ('loop', int)]
    
    @classmethod
    def get_dtype(cls):
        return cls._dty