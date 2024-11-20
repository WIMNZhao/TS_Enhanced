# Enhanced Thompson Sampling for Virtual Screening

This repo accompanies the preprint ["Enhanced Thompson Sampling by Roulette Wheel Selection for
Screening Ultra-Large Combinatorial Libraries"](https://www.biorxiv.org/content/10.1101/2024.05.16.594622v1).

### Setting up the environment for running Enhanced Thompson Sampling

Create a new conda environment and install rdkit:
`conda create -c conda-forge -n <your-env-name> rdkit`

Activate your environment and install the rest of the requirements:
`conda activate <your-env-name>`
`pip install -r requirements.txt`

### Reproduce Figure 3 in the preprint

`python ./src/ts_paper.py quinazoline_fp_sim.json`

Code in src is outdated and will not be maintained

### Use multiprocessing from src_multiprocess

`python ./src_multiprocess/ts_main.py input.json`

Note that there is multiprocessing overhead. For a very efficient scoring method such as fingerprints similarity, it took 1 minute using 1 process while it took 18 minutes using 4 processes (num_per_cylce = 100).

However, by adding time.sleep(0.1) to the evaluation routine, the compuational time scaled linearly with the number of processes. Therefore, for time-limiting scoring such as docking, multiprocess is preferred.
Make a test to check if there is any benefit from multiprocessing before using all CPUs on a cluster, e.g., by setting nprocesses to 1 and 4, respectively, for a comparison.

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
- `num_warmup_trials`: Optional. Default 5. Number of times to randomly sample each reagent in the reaction. It has a big impact on how quickly top-scored compounds can be recovered. Try 5/10/20.
- `results_filename`: Optional. Name of the file to output results to. If None, results will not be saved to a file.
- `log_filename`: Optional. Log filename to save logs to. If not set, logging will be printed to stdout.
- `scaling`: Optional. Default 1. Positive if higher score is preferred; negative otherwise. +/- 1 is usually good. It is used to scale the Boltzmann temperature.
- `decay`: Optional. Default 1. Temperature control
- `num_per_cycle`: Optional. Default 100. A high number helps reduce the TS overhead a bit.
- `stop`: Optional. Default 1000. Stop searching if a new compound has not been sampled for a specified number of consecutive attempts. Increasing the num_warmup_trials may lead to an unexpected early stop. If so, increase it to a higher number.

