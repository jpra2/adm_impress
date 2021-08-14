from packs.data_class import DataManager

class CompositionalData(DataManager):
    _dty = [('prod_o', float), ('prod_w', float), ('wor', float),
       ('delta_t', float), ('dvpi', float), ('loop', int), ('global_identifier', int)]
    
    @classmethod
    def get_dtype(cls):
        return cls._dty
