class MetaDict:
    __slots__ = ['_data']

    def clear(self):
        return self._data.clear()

    def copy(self):
        return self._data.copy()

    def has_key(self, k):
        return k in self._data

    def update(self, *args, **kwargs):
        return self._data.update(*args, **kwargs)

    def keys(self):
        return self._data.keys()

    def values(self):
        return self._data.values()

    def items(self):
        return self._data.items()

    def pop(self, *args):
        return self._data.pop(*args)

    def insert_data(self, data:dict):
        if isinstance(data, dict):
            pass
        else:
            raise TypeError("the data isn't a dict isntance")
        self._data = data

    def __cmp__(self, dict_):
        return self.__cmp__(self._data, dict_)

    def __contains__(self, item):
        return item in self._data

    def __iter__(self):
        return iter(self._data)


    def __str__(self):
        return str(type(self))

    def __setitem__(self, key, value):
        self._data[key] = value

    def __getitem__(self, key):
        return self._data[key]

    def __delitem__(self, key):
        del self._data[key]

    def __hash__(self):
        return hash(self._data)

    def __del__(self):
        del self._data