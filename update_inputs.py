

def update_inputs():
    from packs.utils.info_manager import InfoManager
    import os
    dd = InfoManager('input_cards/inputs0_2.yml', 'input_cards/inputs0.yml')
    dd['load_data'] = True
    dd['load_multilevel_data'] = True
    dd['load_operators'] = True
    dd['load_biphasic_data'] = True

    dd['read_permeability'] = True
    dd['set_permeability'] = False

    # dd['deletar_results'] = True
    dd.save_obj()

    if dd['deletar_results']:

        results = 'results'
        ff = os.listdir(results)

        for f in ff:
            # if f[-4:] == '.vtk':
            #     os.remove(os.path.join(results, f))

            if ff.endswith('.vtk'):
                os.remove(os.path.join(results, f))

update_inputs()
del update_inputs
