from libs import data_loading_maristany, functions
from libs.data_loading_generated import calculate_features_simulation, calculate_features_simulation_v2

def calculate_features_maristany(cfg_path="configs/features_config.json",
                                 data_path="data/seqs_maristany.json"):

    data = data_loading_maristany.load_and_validate_data(data_path)
    cfg = functions.load_cfg(cfg_path)
    return functions.calculate_features(data, cfg)

# data_df = calculate_features_simulation("configs/features_config.json", "data/generated_data")
data_df = calculate_features_simulation("configs/features_config.json", "data/FINISHED_SIMULATIONS", "v3")
print(data_df.shape)


# calculate_features_maristany()

#calculate_features_simulation()

# cfg_path="configs/features_config.json"
# cfg = functions.load_cfg(cfg_path)
# f = functions.calculate_features_from_seqs({"TEST": "ASDSADSADASDASDASD"}, cfg)
# z=1
