import os
import yaml


with open('input_cards/inputs0.yml', 'r') as f:data_loaded = yaml.safe_load(f)
value=float(data_loaded['Permeability']['r2']['value'][0])
import pdb; pdb.set_trace()
while value<10:
    value=10*value
    data_loaded['Permeability']['r2']['value'][0]=value
    data_loaded['Permeability']['r2']['value'][4]=value
    data_loaded['Permeability']['r3']['value'][0]=value
    data_loaded['Permeability']['r3']['value'][4]=value
    with open('input_cards/inputs0.yml', 'w') as file:
        documents = yaml.dump(data_loaded, file, default_flow_style=False)
    os.system('python ADM.py')
