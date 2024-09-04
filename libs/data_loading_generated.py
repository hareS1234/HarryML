<<<<<<< HEAD
import pandas as pd

from libs import functions
=======
from libs import functions
import os
import pandas as pd

>>>>>>> 30d3077 (Extra functions and loaders)


def _generate_df(data_with_features, temp_inc=True):
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
        if temp_inc:
            for temp, c_values in values["temp"].items():
                genes_features[gene_name][str(temp) + "_T"] = c_values['conc']

        for estimator_name, estimator_values in values["features"].items():
            est_abr = est_names_abr[estimator_name]
            for feature_name, feature_value in estimator_values.items():
                genes_features[gene_name][feature_name + f'({est_abr})'] = feature_value

    return pd.DataFrame.from_dict(genes_features, orient='index')


def calculate_features_simulation(cfg_path,
                                  data_path):

    data = functions.load_ps_concentrations(data_path)
    cfg = functions.load_cfg(cfg_path)
    genes_seqs = {gene: values_dict['seq'] for gene, values_dict in data.items()}
    genes_features = functions.calculate_features_from_seqs(genes_seqs, cfg)
    for gene_name in genes_features:
        data[gene_name]["features"] = genes_features[gene_name]["features"]
    return _generate_df(data)

<<<<<<< HEAD
=======

>>>>>>> 30d3077 (Extra functions and loaders)
# data_folder = os.path.join(os.getcwd(), "..\data\generated_data")
# generate_df(data_folder)