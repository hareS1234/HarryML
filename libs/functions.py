import copy
import inspect
import json
import os
from typing import Dict, List, Union

import libs.da as da

import sparrow
from Bio.SeqUtils.ProtParam import ProteinAnalysis


class CustomEstimator:
    def __init__(self, prot_seq: str) -> None:
        self.prot_seq = prot_seq

    def arginine_fraction(self) -> float:
        return sum([1 if b == "R" else 0 for b in self.prot_seq]) / len(self.prot_seq)


estimator_types = Union[CustomEstimator, ProteinAnalysis, sparrow.Protein, sparrow.Protein]
                        

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

def load_ps_concentrations_v2(root_dir):
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

    def parse_critical_temperature(ps_values: list) -> dict:
        for line in ps_values:
            if line.startswith("Critical Solution Temperature"):
                ct_value = line.strip().split(":")[1].strip().split(" ")[0].strip()
        return {"ct_value": float(ct_value)}

    ps_conc = {}
    for gene_name in os.listdir(root_dir):
        ps_conc[gene_name] = {}
        ps_conc[gene_name]["seq"] = _load_seq_from_file(os.path.join(root_dir, gene_name, "seq.txt"))
        ps_conc[gene_name]["temp"] = {}
        for file_name in os.listdir(os.path.join(root_dir, gene_name)):
            if not file_name.startswith(".") and is_chuncked_file_name(file_name):
                temp = file_name.split(".")[0][-3:]
                with open(os.path.join(root_dir, gene_name, file_name)) as f:
                    ps_conc[gene_name]["temp"][float(temp)] = f.readlines()
                    ps_conc[gene_name]["temp"][float(temp)] = parse_ps_concetrations_data(ps_conc[gene_name]["temp"][float(temp)])
            elif not file_name.startswith(".") and file_name == 'critical_solution_temperature.txt':
                with open(os.path.join(root_dir, gene_name, file_name)) as f:
                    ps_conc[gene_name]["ct"] = f.readlines()
                    ps_conc[gene_name]["ct"] = parse_critical_temperature(ps_conc[gene_name]["ct"])
    return ps_conc


def load_ps_concentrations_v3(root_dir):
    def is_chuncked_file_name(file_name):
        if file_name.endswith(".txt") and "densities_chunked2" in file_name:
            return True

        #if file_name.endswith(".txt") and file_name.split(".")[0][-3:].isdigit():
        #    return True
        return False
    
    def is_raw_chuncked_file_name(file_name):
        if file_name.endswith(".txt") and file_name.split(".")[0][-3:].isdigit():
            return True
        return False

    def parse_ps_concetrations_data(ps_values: list) -> dict:
        
        top_key = 'conc'
        d_values = {top_key: {}}
        
        for line in ps_values:
            values = line.strip().split(":")[1].split("+/-")
            conc = float(values[0].strip())
            var = float(values[1].strip())
            
            if line.startswith("Dense"):
               d_values[top_key]['dense_conc'] = conc
               d_values[top_key]['var_dense'] = var
        
            elif line.startswith("Dilute"):
                d_values[top_key]['conc_dilute'] = conc
                d_values[top_key]['var_dilute'] = var
        
        return d_values

    def parse_critical_temperature(ps_values: list) -> dict:
        
        for line in ps_values:
        
            if line.startswith("Critical Solution Temperature"):
                ct_value = line.strip().split(":")[1].strip().split(" ")[0].strip()
        
                if "-" in ct_value:
                    ct_value = ct_value.split("-")[0]
        
        return {"ct_value": float(ct_value)}

    ps_conc = {}
    
    for gene_name in os.listdir(root_dir):
        ps_conc[gene_name] = {}
        ps_conc[gene_name]["seq"] = _load_seq_from_file(os.path.join(root_dir, gene_name, "seq.txt"))
        ps_conc[gene_name]["temp"] = {}
        # ps_conc[gene_name]["density"] = {}
        # ps_conc[gene_name]["stddevs"] = {}

        for temp in os.listdir(os.path.join(root_dir, gene_name)):
            
            if temp[4:].isdigit() or temp.isdigit():
                
                if temp.isdigit():
                    temp_v = float(temp)
                else:
                    temp_v = float(temp[4:])
                
                if float(temp_v) not in ps_conc[gene_name]["temp"]:
                    ps_conc[gene_name]["temp"][temp_v] = {}
                
                for file_name in os.listdir(os.path.join(root_dir, gene_name, temp)):
                
                    if not file_name.startswith(".") and is_chuncked_file_name(file_name):
                        #temp_v = file_name.split(".")[0][-3:]
                        with open(os.path.join(root_dir, gene_name, temp, file_name)) as f:
                            data_t = f.readlines()
                            #ps_conc[gene_name]["temp"][float(temp_v)] = f.readlines()
                            ps_conc[gene_name]["temp"][float(temp_v)]["conc"] = parse_ps_concetrations_data(data_t)["conc"]

                    elif not file_name.startswith(".") and file_name == 'densities_chunked2.dat':
                        
                        density, stddevs = da.get_density(os.path.join(root_dir, gene_name, temp, file_name))
                        ps_conc[gene_name]["temp"][float(temp_v)]["density"] = list(density)
                        ps_conc[gene_name]["temp"][float(temp_v)]["stddevs"] = list(stddevs)
                        #ps_conc[gene_name]["temp"][float(temp_v)] = stddevs
            
            elif not temp.startswith(".") and temp == 'critical_solution_temperature.txt':
            
                with open(os.path.join(root_dir, gene_name, temp)) as f:
                    ps_conc[gene_name]["ct"] = f.readlines()
                    ps_conc[gene_name]["ct"] = parse_critical_temperature(ps_conc[gene_name]["ct"])
            
            elif not temp.startswith(".") and temp == 'ps_densities':
            
                for file_to_check in os.listdir(os.path.join(root_dir, gene_name, temp)):
            
                    if file_to_check == 'critical_solution_temperature.txt':
            
                        with open(os.path.join(root_dir, gene_name, temp, file_to_check)) as f:
                            ps_conc[gene_name]["ct"] = f.readlines()
                            ps_conc[gene_name]["ct"] = parse_critical_temperature(ps_conc[gene_name]["ct"])
    return ps_conc


# TODO add this
def sequence_identity_and_similarity():
    pass


def _calculate_features(
    estimator: estimator_types, 
    features_to_calculate: Dict[str, list]
) -> Dict[str, float]:
    """"""
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


def _calculate_features_all(prot_seq: str, cfg: dict, ESTIMATORS: Dict[str, estimator_types]):
    # TODO change the way predictor put in config
    features = {}
    
    for estimator_name, estimator_features in cfg["ESTIMATORS"].items():
        
        if estimator_name != "SPARROW_PREDICTOR":
            estimator = ESTIMATORS[estimator_name](prot_seq)
        else:
            estimator = ESTIMATORS[estimator_name](prot_seq).predictor
            
        features[estimator_name] = _calculate_features(
            estimator=estimator,
            features_to_calculate=estimator_features
        )
        
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


def calculate_features_from_seqs(seqs: List[str], cfg: dict):
    
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

















