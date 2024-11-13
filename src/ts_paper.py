#!/usr/bin/env python


import sys
from datetime import timedelta
from timeit import default_timer as timer

from ts_utils import read_input
from run_ts import run_ts


def run_10_replicates():
    """ 
    A testing function for the paper fig.4
    :return: None
    """
    json_file_name = sys.argv[1]
    input_dict = read_input(json_file_name)
    for i in range(0, 10):
        input_dict['results_filename'] = f"./results/ts_result_1_{i:03d}.csv"
        input_dict['num_per_cycle'] = 1
        input_dict['num_ts_iterations'] = 210000
        run_ts(input_dict, hide_progress=False)

    for i in range(0, 10):
        input_dict['results_filename'] = f"./results/ts_result_3_{i:03d}.csv"
        input_dict['num_per_cycle'] = 3
        input_dict['num_ts_iterations'] = 70000
        run_ts(input_dict, hide_progress=False)


if __name__ == "__main__":
    run_10_replicates()
