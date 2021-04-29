def test_kwargs_keys(default_keys, keys):
    resp = set(keys) & default_keys
    
    if resp == keys:
        pass
    else:
        raise ValueError(f'These kwargs keys not in default keys.\n \
                         Default kwargs keys: {default_keys}.\n \
                         kwargs keys: {keys} \n')