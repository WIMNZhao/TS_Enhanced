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

Note that there is multiprocessing overhead. If time_for_scoring_single_compound * num_per_cycle / nprocessses > 0.1, there is likely a gain. For example, it took 34 seconds to screening 0.1% the 94M quinazoline library using 1 process (num_per_cycle = 1000). While it took 26 seconds using 2 processes. For expensive scoring (e.g., by add time.sleep(0.1) to the evaluaiton routine), it will scale linearly withe CPU cores.

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
- `num_warmup_trials`: Optional. Default 5. Number of times to randomly sample each reagent in the reaction. It has an impact on how quickly top-scored compounds can be recovered. For the size of the largest reaction component library around a few hundreds, try 20. For a size around tens of thousands, try 5.
- `num_per_cycle`: Optional. Default 100. size_of_largest_component_library * n with n in the range from 1 to 4. For instance, for the quinazoline library, the search with 2000 compounds per cycle significantly outperforms that with 2 compounds per cycyle. A high number could also benefit from multiprocessing a lot. 
- `scaling`: Optional. Default 1. Positive if higher score is preferred; negative otherwise. +/- 1 is usually good. It is used to scale the Boltzmann temperature.
- `decay`: Optional. Default 1. Temperature control
- `stop`: Optional. Default 1000. Stop searching when without sampling a new compound for a specified number of consecutive attempts. Increasing the num_warmup_trials may lead to an unexpected early stop. If so, increase it to a higher number.
- `results_filename`: Optional. Name of the file to output results to. If None, results will not be saved to a file.
- `log_filename`: Optional. Log filename to save logs to. If not set, logging will be printed to stdout.

