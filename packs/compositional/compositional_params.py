from packs.data_class.meta_dict import MetaDict

class Params(MetaDict):
    def __init__(self):
        self.insert_data(dict())
    pass