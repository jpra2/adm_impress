
class CommonData1:

    def get_value(self, name):
        return self.data[name]

    def set_value(self, name, data):
        self.data[name] = data

    def datas_info(self):
        return self.data.keys()
