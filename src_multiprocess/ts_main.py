#!/usr/bin/env python

import sys
from datetime import timedelta
from timeit import default_timer as timer
from ts_utils import read_input
from ts_run import run_ts
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')


def main():
    start = timer()
    json_filename = sys.argv[1]
    input_dict = read_input(json_filename)
    run_ts(input_dict)
    end = timer()
    print("Elapsed time", timedelta(seconds=end - start))


if __name__ == "__main__":
    main()
