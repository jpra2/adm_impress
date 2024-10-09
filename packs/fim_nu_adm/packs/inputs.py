import os
import yaml
parent_dir=os.getcwd()
with open(parent_dir+'/inputs/finescale.yml', 'r') as f:
    finescale_inputs = yaml.safe_load(f)

with open(parent_dir+'/inputs/multiscale_and_multilevel.yml', 'r') as f:
    multiscale_and_multilevel_inputs = yaml.safe_load(f)
