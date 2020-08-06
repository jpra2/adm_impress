
class CommonDataManager:

    def keys(self):
        return self._data.keys()

    def values(self):
        return self._data.values()

    def items(self):
        return self._data.items()

    def removekey(self, d, key):
        r = dict(d)
        del r[key]
        return r

    def update(self, data):
        self._data.update(data)

    def __setitem__(self, key, value):
        self._data[key] = value

    def __getitem__(self, key):
        return self._data[key]

    def __delitem__(self, key):
        del self._data[key]

    def __hash__(self, key):
        return hash(self._data)

    def __contains__(self, key):
        return key in self._data
