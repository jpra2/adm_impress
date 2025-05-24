import os

def run():
    

    list_packages = [
        'pandas==2.0.3',
        'h5sparse==0.1.0',
        'meshio==5.3.5',
        'shapely==2.0.7',
        'Pint==0.21.1',
        'matplotlib==3.7.5'
    ]

    for name in list_packages:
        instruction = 'pip install ' + name
        os.system(instruction)