from xml.dom import ValidationErr
from packs.data_class.meta_dict import MetaDict

class Params(MetaDict):
    def __init__(self):
        self.insert_data(dict())
    
    def get(self, key):
        return self._data[key]