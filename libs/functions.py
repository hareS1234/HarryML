import copy
import inspect
import json
import os

import sparrow
from Bio.SeqUtils.ProtParam import ProteinAnalysis


def _get_prot_biopython(prot_seq: str) -> ProteinAnalysis:
    return ProteinAnalysis(prot_seq)


def _get_prot_sparrow(prot_seq: str) -> sparrow.Protein:
    return sparrow.Protein(prot_seq)


def _get_prot_sparrow_pred(prot_seq: str) -> sparrow.Protein.predictor:
    return sparrow.Protein(prot_seq)


def load_cfg(file_path):
    with open(file_path) as f:
        return json.load(f)


def _load_seq_from_file(file_name="seq.txt"):
    with open(file_name) as f:
        seq = f.readlines()
    return "".join([l.strip() for l in seq])


def load_ps_concentrations(root_dir):
    def is_chuncked_file_name(file_name):
        if file_name.endswith(".txt") and file_name.split(".")[0][-3:].isdigit():
            return True
        return False

    def parse_ps_concetrations_data(ps_values: list) -> dict:
        for line in ps_values:
            if line.startswith("Dense"):
                values = line.strip().split(":")[1].split("+/-")
                conc_dense = float(values[0].strip())
                var_dense = float(values[1].strip())
            elif line.startswith("Dilute"):
                values = line.strip().split(":")[1].split("+/-")
                conc_dilute = float(values[0].strip())
                var_dilute = float(values[1].strip())

        return {"conc": {"dilute_conc": conc_dilute,
                         "dense_conc": conc_dense,
                         "var_dilute": var_dilute,
                         "var_dense": var_dense}}

    ps_conc = {}
    for gene_name in os.listdir(root_dir):
        ps_conc[gene_name] = {}
        ps_conc[gene_name]["seq"] = _load_seq_from_file(os.path.join(root_dir, gene_name, "seq.txt"))
        ps_conc[gene_name]["temp"] = {}
        for temp in os.listdir(os.path.join(root_dir, gene_name)):
            if temp.split(os.sep)[-1].isdigit():
                for file_name in os.listdir(os.path.join(root_dir, gene_name, temp)):
                    if is_chuncked_file_name(file_name):
                        with open(os.path.join(root_dir, gene_name, temp, file_name)) as f:
                            ps_conc[gene_name]["temp"][float(temp)] = f.readlines()
                            ps_conc[gene_name]["temp"][float(temp)] = parse_ps_concetrations_data(ps_conc[gene_name]["temp"][float(temp)])
    return ps_conc


class CustomEstimator:
    def __init__(self, prot_seq: str) -> None:
        self.prot_seq = prot_seq

    def arginine_fraction(self) -> float:
        return sum([1 if b == "R" else 0 for b in self.prot_seq]) / len(self.prot_seq)


# TODO add this
def sequence_identity_and_similarity():
    pass


def _calculate_features(estimator, features_to_calculate):
    if len(list(features_to_calculate.keys())) != len(set(list(features_to_calculate.keys()))):
        raise ValueError(f"Duplicate keys found for {estimator}")

    features_calculated = {}
    for feature_name, (function_name, args, kwargs) in features_to_calculate.items():
        func = getattr(estimator, function_name)
        if inspect.ismethod(func):
            features_calculated[feature_name] = func(*args, **kwargs)
        else:
            if args:
                raise ValueError(f"non-empty arguments for class property - {function_name}")
            else:
                features_calculated[feature_name] = func

    return features_calculated

def _calculate_features_all(prot_seq, cfg, ESTIMATORS):
    # TODO change the way predictor put in config
    features = {}
    for estimator_name, estimator_features in cfg["ESTIMATORS"].items():
        if estimator_name != "SPARROW_PREDICTOR":
            estimator = ESTIMATORS[estimator_name](prot_seq)
        else:
            estimator = ESTIMATORS[estimator_name](prot_seq).predictor
        features[estimator_name] = _calculate_features(estimator=estimator,
                                                       features_to_calculate=estimator_features)
    return features

def _calculate_features_all(prot_seq, cfg, ESTIMATORS):
    # TODO change the way predictor put in config
    features = {}
    for estimator_name, estimator_features in cfg["ESTIMATORS"].items():
        if estimator_name != "SPARROW_PREDICTOR":
            estimator = ESTIMATORS[estimator_name](prot_seq)
        else:
            estimator = ESTIMATORS[estimator_name](prot_seq).predictor
        features[estimator_name] = _calculate_features(estimator=estimator,
                                                       features_to_calculate=estimator_features)
    return features


def calculate_features(data, cfg):
    data_with_features = copy.deepcopy(data)
    ESTIMATORS = {
        "CUSTOM": CustomEstimator,
        "BIOPYTHON": ProteinAnalysis,
        "SPARROW": sparrow.Protein,
        "SPARROW_PREDICTOR": sparrow.Protein
    }

    for i, entry in enumerate(data_with_features):
        print(f"processing block {i} - gene: {entry['gene']}")
        entry["seq_base"] = {
            "seq_base": entry["seq_base"],
            "features": _calculate_features_all(entry["seq_base"], cfg, ESTIMATORS)
        }
        mut_data_with_features = {}
        for mut_name, mut_seq in entry["seq_gen"].items():
            mut_data_with_features[mut_name] = {
                "mut_base": mut_seq,
                "features": _calculate_features_all(mut_seq)
            }
        entry["seq_gen"] = mut_data_with_features

    return data_with_features


def calculate_features_from_seqs(seqs, cfg):
    ESTIMATORS = {
        "CUSTOM": CustomEstimator,
        "BIOPYTHON": ProteinAnalysis,
        "SPARROW": sparrow.Protein,
        "SPARROW_PREDICTOR": sparrow.Protein
    }
    seqs_with_features = {}
    for gene, seq in seqs.items():
        seqs_with_features[gene] = {}
        seqs_with_features[gene]["seq"] = seq
        seqs_with_features[gene]["features"] = _calculate_features_all(seq, cfg, ESTIMATORS)
    return seqs_with_features

















