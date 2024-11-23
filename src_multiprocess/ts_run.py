import pandas as pd
from thompson_sampling import ThompsonSampler
from ts_logger import get_logger


def run_ts(input_dict: dict) -> None:
    """
    Run a Thompson sampling with roulette wheel selection
    :param hide_progress: hide the progress bar
    :param input_dict: dictionary with input parameters
    """
    search = {
        "percent_of_library": None,
        "num_per_cycle": None,
        "scaling": None,
        "stop": None,
        "results_filename": None
    }
    for para in search:
        search[para] = input_dict[para]

    logger = get_logger(__name__, filename=input_dict["log_filename"])

    # setup ts
    ts = ThompsonSampler(input_dict["nprocesses"])
    ts.set_hide_progress(input_dict["hide_progress"])
    ts.set_evaluator(input_dict["evaluator_class"])
    ts.read_reagents(input_dict["reagent_file_list"], num_to_select=None)
    ts.set_reaction(input_dict["reaction_smarts"])
    # run the warm-up phase to generate an initial set of scores for each reagent
    nw = ts.warm_up(input_dict["num_warmup_trials"], input_dict["results_filename"])
    # run the search with TS
    out_list = ts.search(**search)

    # logging
    total_evaluations = len(out_list) + nw
    percent_searched = total_evaluations / ts.get_num_prods() * 100
    logger.info(f"{total_evaluations} evaluations | {percent_searched:.3f}% of total")
    logger.info(f"Saved results to: " + input_dict["results_filename"])
    out_df = pd.DataFrame(out_list, columns=["score", "SMILES", "Name"])
    if not input_dict["hide_progress"]:
        if input_dict["scaling"] > 0:
           print(out_df.sort_values("score", ascending=False).drop_duplicates(subset="SMILES").head(100))    
        else:
           print(out_df.sort_values("score", ascending=True).drop_duplicates(subset="SMILES").head(100))  
    
