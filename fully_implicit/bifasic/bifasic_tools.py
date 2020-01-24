from pyaml import yaml

def set_wells():
    with open("input/wells.yml", 'r') as stream:
        wells_loaded = yaml.safe_load(stream)
