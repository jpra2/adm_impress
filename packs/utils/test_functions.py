def test_kwargs_keys(default_keys, keys):
    
    # resp = set(keys) & set(default_keys)
    
    # if resp == keys:
    if set(keys) & set(default_keys) == keys:
        pass
    else:
        raise ValueError(f'These kwargs keys not in default keys.\n \
                         Default kwargs keys: {default_keys}.\n \
                         kwargs keys: {keys} \n')


def test_instance(obj, inst):
    if isinstance(obj, inst):
        pass
    else:
        raise TypeError(f'The object {obj} not is a instance of {inst}')