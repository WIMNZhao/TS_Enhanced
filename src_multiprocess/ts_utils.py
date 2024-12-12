from typing import List, Optional
import json, importlib
import multiprocessing as mp
from reagent import Reagent


def create_reagents(filename: str, num_to_select: Optional[int] = None) -> List[Reagent]:
    """
    Creates a list of Reagents from a file
    :param filename: a smiles file containing the reagents
    :param num_to_select: For dev purposes; the number of molecules to return
    :return: List of Reagents
    """
    reagent_list = []
    with open(filename, 'r') as f:
        for line in f.readlines():
            smiles, reagent_name = line.split()
            reagent = Reagent(reagent_name=reagent_name, smiles=smiles)
            reagent_list.append(reagent)
    if num_to_select is not None and len(reagent_list) > num_to_select:
        reagent_list = reagent_list[:num_to_select]
    return reagent_list

def read_reagents(reagent_file_list, num_to_select: Optional[int]) -> List[Reagent]:
    """
    Read the reagents SMILES files
    :param reagent_file_list: a list of filenames containing reagents for the reaction. Each file list contains smiles
                              strings for a single component of the reaction.
    :param num_to_select: select how many reagents to read, mostly a development function
    :return: List of Reagents
    """
    reagents = []
    for reagent_filename in reagent_file_list:
        reagent_list = create_reagents(filename=reagent_filename, num_to_select=num_to_select)
        reagents.append(reagent_list)
    return reagents

def read_input(json_filename: str) -> dict:
    """
    Read input parameters from a json file
    :param json_filename: input json file
    :return: a dictionary with the input parameters
    """
    input_data = None
    with open(json_filename, 'r') as ifs:
        input_data = json.load(ifs)
        module = importlib.import_module("evaluators")
        evaluator_class_name = input_data["evaluator_class_name"]
        class_ = getattr(module, evaluator_class_name)
        evaluator_arg = input_data["evaluator_arg"]
        evaluator = class_(evaluator_arg)
        input_data['evaluator_class'] = evaluator
    default = {
        "nprocesses": mp.cpu_count(),
        "num_warmup_trials": 3,
        "percent_of_library": 0.001,
        "num_per_cycle": 1000,
        "scaling": 1,
        "stop": 1000,
        "hide_progress": True,
        "log_filename": "./logs.txt",
        "results_filename": "./results.csv"
    }
    for para in default:
        if para not in input_data:
           input_data[para] = default[para]
    return input_data

