
def get_config_preset_list_for_model(weight_set, model_name):
    if weight_set == "OpenFold" and model_name == "monomer":
        model_names = [
            'finetuning_3.pt',
            'finetuning_4.pt',
            'finetuning_5.pt',
            'finetuning_ptm_2.pt',
            'finetuning_no_templ_ptm_1.pt'
        ]
    elif weight_set == "AlphaFold" and model_name == "multimer":
        # Generate 'model_1_multimer_v3' to 'model_5_multimer_v3' using a loop
        model_names = [f'model_{i}_multimer_v3' for i in range(1, 6)]
    elif weight_set == "AlphaFold" and model_name == "monomer":
        # Generate 'model_1' to 'model_5' using a loop
        model_names = [f'model_{i}' for i in range(1, 6)]
    else:
        raise ValueError(f"Invalid combination of weight_set '{weight_set}' and model_name '{model_name}'")
    
    return model_names