import pandas as pd
import numpy as np
from libs import data_loading_maristany, functions, data_loading_generated
from libs.data_loading_generated import calculate_features_simulation, calculate_features_simulation_v2
import matplotlib.pyplot as plt

from sklearn.svm import SVR
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split

from sklearn.tree import DecisionTreeRegressor
import copy

from sklearn.model_selection import GridSearchCV

from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import Normalizer

from sklearn.metrics import mean_squared_error

import math

def clear_df(data_df):
    columns_to_drop = []
    for column_name in data_df.columns:
        if not isinstance(data_df.iloc[0][column_name], np.number) and column_name != "seq" and  not column_name.endswith("_T"):
            columns_to_drop.append(column_name)
    data_df = data_df.drop(labels=columns_to_drop, axis=1)

data_df = calculate_features_simulation("configs/features_config.json", "data/FINISHED_SIMULATIONS/", calc_type="v3",
                                        dens_raw=False)
columns_to_drop = []
for column_name in data_df.columns:
    if not isinstance(data_df.iloc[0][column_name], np.number) and column_name != "seq" and  not column_name.endswith("_T"):
        columns_to_drop.append(column_name)
data_df = data_df.drop(labels=columns_to_drop, axis=1)