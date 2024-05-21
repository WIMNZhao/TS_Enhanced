# Enhanced Thompson Sampling for Virtual Screening

This repo accompanies the preprint ["Enhanced Thompson Sampling by Roulette Wheel Selection for
Screening Ultra-Large Combinatorial Libraries"](https://www.biorxiv.org/content/10.1101/2024.05.16.594622v1).

### Setting up the environment for running Enhanced Thompson Sampling

Create a new conda environment and install rdkit:
`conda create -c conda-forge -n <your-env-name> rdkit`

Activate your environment and install the rest of the requirements:
`conda activate <your-env-name>`
`pip install -r requirements.txt`

### How to run Thompson Sampling

`python ./src/ts_main.py test.json`

Or try the example used in the preprint:

`python ./src/ts_paper.py quinazoline_fp_sim.json`

### Parameters

Required params:
- `evaluator_arg`: The argument to be passed to the instantiation of the evaluator (e.g. smiles string of the query
molecule, filepath to the mol file, etc.) See the relevant `Evaluator` for more info.
- `evaluator_class_name`: Required. The scoring function to use. Must be "FPEvaluator" for 2D similarity. To use your own scoring function, implement a subclass of the
Evaluator baseclass in evaluators.py.
- `reaction_smarts`: Required. The SMARTS string for the reaction.
- `num_ts_iterations`: Required. Number of iterations of Thompson Sampling to run (usually 100 - 2000).
- `reagent_file_list`: Required. List of filepaths - one for each component of the reaction. Each file should contain the
smiles strings of valid reagents for that component.

Optional params:
- `num_warmup_trials`: Optional. Number of times to randomly sample each reagent in the reaction. 3 is usually sufficient
- `results_filename`: Optional. Name of the file to output results to. If None, results will not be saved to a file.
- `log_filename`: Optional. Log filename to save logs to. If not set, logging will be printed to stdout.
- `scaling`: Optional. Positive if higher score is preferred; negative otherwise. +/- 1 is usually good
- `decay`: Optional. Temperature control
- `num_per_cycle`: Optional. e.g., 3, 10, 100

