import pandas as pd

from thompson_sampling import ThompsonSampler
from ts_logger import get_logger


def run_ts(input_dict: dict, hide_progress: bool = False) -> None:
    """
    Run a Thompson sampling with roulette wheel selection
    :param hide_progress: hide the progress bar
    :param input_dict: dictionary with input parameters
    """
    evaluator = input_dict["evaluator_class"]
    reaction_smarts = input_dict["reaction_smarts"]
    num_ts_iterations = input_dict["num_ts_iterations"]
    reagent_file_list = input_dict["reagent_file_list"]
    num_warmup_trials = input_dict["num_warmup_trials"]
    scaling = input_dict["scaling"]
    decay = input_dict["decay"]
    num_per_cycle = input_dict["num_per_cycle"]
    result_filename = input_dict.get("results_filename")
    log_filename = input_dict.get("log_filename")
    logger = get_logger(__name__, filename=log_filename)

    ts = ThompsonSampler()
    ts.set_hide_progress(hide_progress)
    ts.set_evaluator(evaluator)
    ts.read_reagents(reagent_file_list=reagent_file_list, num_to_select=None)
    ts.set_reaction(reaction_smarts)
    # run the warm-up phase to generate an initial set of scores for each reagent
    ts.warm_up(num_warmup_trials=num_warmup_trials)
    # run the search with TS
    out_list = ts.search(scaling=scaling,decay=decay,num_per_cycle=num_per_cycle,num_cycles=num_ts_iterations)

    total_evaluations = evaluator.counter
    percent_searched = total_evaluations / ts.get_num_prods() * 100
    logger.info(f"{total_evaluations} evaluations | {percent_searched:.3f}% of total")
    # write the results to disk
    out_df = pd.DataFrame(out_list, columns=["score", "SMILES", "Name"])
    if result_filename is not None:
        out_df.to_csv(result_filename, index=False)
        logger.info(f"Saved results to: {result_filename}")
    if not hide_progress:
        if scaling > 0:
           print(out_df.sort_values("score", ascending=False).drop_duplicates(subset="SMILES").head(100))    
        else:
           print(out_df.sort_values("score", ascending=True).drop_duplicates(subset="SMILES").head(100))  
    