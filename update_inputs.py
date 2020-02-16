

def update_inputs():
    from packs.utils.info_manager import InfoManager
    import os
    dd = InfoManager('input_cards/inputs0_2.yml', 'input_cards/inputs0.yml')
    # dd['load_data'] = True
    # dd['load_multilevel_data'] = True
    # dd['load_operators'] = True
    # dd['load_biphasic_data'] = True

    # dd['read_permeability'] = True
    # dd['set_permeability'] = False

    dd['biphasic'] = False

    dd['deletar_results'] = True
    # dd['gravity'] = True
    # dd['get_correction_term'] = True
    # dd['convert_english_to_SI'] = True
    dd.save_obj()

    dd2 = InfoManager('input_cards/variable_adm.yml', 'input_cards/variable_input.yml')
    dd2.save_obj()

    if dd['deletar_results']:

        results = 'results'
        ff = os.listdir(results)

        for f in ff:
            # if f[-4:] == '.vtk':
            #     os.remove(os.path.join(results, f))

            if f.endswith('.vtk'):
                os.remove(os.path.join(results, f))

update_inputs()
del update_inputs
