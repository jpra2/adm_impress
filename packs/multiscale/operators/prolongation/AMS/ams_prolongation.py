from .....data_class.data_manager import DataManager

class AMSTpfa(DataManager):
    name = 'AMSTpfa_'
    id = 1

    def __init__(self, T: 'transmissibility_matrix',
        gids: 'global_ids',
        reordered_ids: 'ids_reordenados',
        load=False):
        
        data_name = AMSTpfa.name + str(AMSTpfa.id)
        super().__init__(data_name=data_name, load=load)
        AMSTpfa.id += 1
