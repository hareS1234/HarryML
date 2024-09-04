import pandas as pd

from libs import functions


def _generate_df(data_with_features, temp_inc=True, ct=True, dens_raw=True):
    est_names_abr = {
        "CUSTOM": "C",
        "BIOPYTHON": "BP",
        "SPARROW": "S",
        "SPARROW_PREDICTOR": "SP"
    }
    genes_features = {}
    for gene_name, values in data_with_features.items():
        genes_features[gene_name] = {}
        genes_features[gene_name]["seq"] = values["seq"]
        # if temp_inc or dens_raw:
        #     genes_features[gene_name][str(temp) + "_T"] = {}
        if temp_inc:
            for temp, c_values in values["temp"].items():
                if str(temp) + "_T" not in genes_features[gene_name]:
                    genes_features[gene_name][str(temp) + "_T"] = {}
                genes_features[gene_name][str(temp) + "_T"]["conc"] = c_values['conc']

        if ct:
            genes_features[gene_name]['ct'] = values['ct']['ct_value']

        if dens_raw:
            for temp, c_values in values["temp"].items():
                if str(temp) + "_T" not in genes_features[gene_name]:
                    genes_features[gene_name][str(temp) + "_T"] = {}
                #if 'density' in c_values:
                genes_features[gene_name][str(temp) + "_T"]["density"] = c_values['density']
                genes_features[gene_name][str(temp) + "_T"]["stddevs"] = c_values['stddevs']

        for estimator_name, estimator_values in values["features"].items():
            est_abr = est_names_abr[estimator_name]
            for feature_name, feature_value in estimator_values.items():
                genes_features[gene_name][feature_name + f'({est_abr})'] = feature_value

    return pd.DataFrame.from_dict(genes_features, orient='index')


def calculate_features_simulation(cfg_path,
                                  data_path,
                                  calc_type="v1",
                                  dens_raw=True):
    if calc_type == 'v1':
        data = functions.load_ps_concentrations(data_path)
    elif calc_type == 'v2':
        data = functions.load_ps_concentrations_v2(data_path)
    elif calc_type == 'v3':
        data = functions.load_ps_concentrations_v3(data_path)
    cfg = functions.load_cfg(cfg_path)
    genes_seqs = {gene: values_dict['seq'] for gene, values_dict in data.items()}
    genes_features = functions.calculate_features_from_seqs(genes_seqs, cfg)
    for gene_name in genes_features:
        data[gene_name]["features"] = genes_features[gene_name]["features"]
    return _generate_df(data, dens_raw=dens_raw)

def calculate_features_simulation_v2(cfg_path,
                                    data_path):

    data = functions.load_ps_concentrations_v2(data_path)
    cfg = functions.load_cfg(cfg_path)
    genes_seqs = {gene: values_dict['seq'] for gene, values_dict in data.items()}
    genes_features = functions.calculate_features_from_seqs(genes_seqs, cfg)
    for gene_name in genes_features:
        data[gene_name]["features"] = genes_features[gene_name]["features"]
    return _generate_df(data)

#data_folder = os.path.join(os.getcwd(), "..\data\generated_data")
#generate_df(data_folder)