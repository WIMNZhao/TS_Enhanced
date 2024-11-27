# Enhanced Thompson Sampling for Virtual Screening

This repo accompanies the preprint ["Enhanced Thompson Sampling by Roulette Wheel Selection for
Screening Ultra-Large Combinatorial Libraries"](https://www.biorxiv.org/content/10.1101/2024.05.16.594622v1).

### Setting up the environment for running Enhanced Thompson Sampling

Create a new conda environment:
`conda create -n ts_test python=3.11`

Activate your environment and install the rest of the requirements:
`conda activate ts_test`
`pip install -r requirements.txt`

### Reproduce Figure 3 in the preprint

`python ./src/ts_paper.py quinazoline_fp_sim.json`

Code in src is outdated and will not be maintained

### Use multiprocessing from src_multiprocess

`python ./src_multiprocess/ts_main.py input.json`

Note that there is a multiprocessing overhead. If time_for_scoring_single_compound * num_per_cycle / nprocessses > e.g., 0.1 s, there is a gain from multiprocessing. It scales linearly for docking!

### Parameters

Required params:
- `evaluator_arg`: The argument to be passed to the instantiation of the evaluator (e.g. smiles string of the query
molecule, filepath to the mol file, etc.) See the relevant `Evaluator` for more info.
- `evaluator_class_name`: Required. The scoring function to use. Must be "FPEvaluator" for 2D similarity. To use your own scoring function, implement a subclass of the
Evaluator baseclass in evaluators.py.
- `reaction_smarts`: Required. The SMARTS string for the reaction.
- `reagent_file_list`: Required. List of filepaths - one for each component of the reaction. Each file should contain the
smiles strings of valid reagents for that component.

Optional params:
- `percent_of_libray`: Optional. Default 0.1%. Percent of library to be screened. If 0.1% to be screened, set as 0.001 instead of 0.1 (a bit confusing).
- `num_warmup_trials`: Optional. Default 5.
- `num_per_cycle`: Optional. Default 1000. A value that is around 10% of the size of the largest component library.
- `scaling`: Optional. Default 1. +1 if higher score is preferred; -1 otherwise. 
- `stop`: Optional. Default 1000. Stop searching when without sampling a new compound for a specified number of consecutive attempts.
- `results_filename`: Optional. Default "./results.csv". Name of the file to output results to.
- `log_filename`: Optional. Log filename to save logs to. If not set, logging will be printed to stdout.
- `hide_progress`: Optional. Defaut true. false otherwise.
- "nprocesses": Optional. Default to use all CPU cores. 

