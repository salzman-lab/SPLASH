#!/usr/bin/env python3
# This script is based on the example at: https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv

import os
from os.path import exists
import sys
import errno
import argparse


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--samplesheet",
        type=str
    )
    parser.add_argument(
        "--is_10X",
        action='store_true'
    )
    args = parser.parse_args()
    return args


def print_error(error, context="Line", context_str=""):
    error_str = "ERROR: Please check samplesheet -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please check samplesheet -> {}\n{}: '{}'".format(
            error, context.strip(), context_str.strip()
        )
    print(error_str)
    sys.exit(1)


def main():
    args = get_args()

    with open(args.samplesheet, "r") as file_in:
        for line in file_in:

            ## get row
            row = [x for x in line.strip().split(',')]

            ## check that samplesheet is formatted as: sample_ID, fastq_R1, fastq_R2
            if args.is_10X:

                ## error if there are not exactly 3 columns
                if len(row) != 3:
                    print_error(
                        "ERROR: 10X samplesheets must contain 3 columns.",
                        "Line",
                        line
                    )

                sample_ID, fastq_R1, fastq_R2 = row

                ## error if any value is empty
                if not sample_ID or not fastq_R1 or not fastq_R2:
                    print_error(
                        f"Error: Column entries in samplesheet must not be empty, for values {row}.",
                        "Line",
                        line
                    )

                ## error if fastqs are not gzipped or if the paths are not valid
                for fastq in [fastq_R1, fastq_R2]:
                    if not fastq.endswith(".fastq.gz") and not fastq.endswith(".fq.gz"):
                        print_error(
                            f"ERROR: FASTQ file does not have extension '.fastq.gz' or '.fq.gz': {fastq}",
                            "Line",
                            line,
                        )
                    if not exists(fastq):
                        print_error(
                            f"ERROR: FASTQ file does not exist at this path: {fastq}",
                            "Line",
                            line,
                        )


            ## check that samplesheet is formatted as: fastq, optional group_ID(-1/1)
            else:
                ## error if there are not exactly 3 columns
                if len(row) > 2:
                    print_error(
                        "ERROR: Samplesheets must contain at most 2 columns.",
                        "Line",
                        line
                    )
                else:
                    ## if no group_ids are provided
                    if len(row) == 1:
                        (fastq_file) = row

                        ## error if any value is empty
                        if not fastq_file:
                            print_error(
                                f"Error: Column entries in samplesheet must not be empty, for values {row}.",
                                "Line",
                                line
                            )
                    ## if group_ids are provided
                    elif len(row) == 2:
                        (fastq_file, group_id) = row

                        ## error if any value is empty
                        if not fastq_file or not group_id:
                            print_error(
                                f"Error: Column entries in samplesheet must not be empty, for values {row}.",
                                "Line",
                                line
                            )

                        ## error if group_id is not 1 or -1
                        if group_id.strip() not in ["-1", "1"]:
                            print_error(
                                f"ERROR: Group ID values must be 1 or -1: {group_id}",
                                "Line",
                                line,
                            )

                    ## error if fastqs are not gzipped or if the paths are not valid
                    try:
                        if not fastq_file.endswith(".fastq.gz") and not fastq_file.endswith(".fq.gz"):
                            print_error(
                                f"ERROR: FASTQ file does not have extension '.fastq.gz' or '.fq.gz': {fastq_file}",
                                "Line",
                                line,
                            )
                        if not exists(fastq_file):
                            print_error(
                                f"ERROR: FASTQ file does not exist at this path: {fastq_file}",
                                "Line",
                                line,
                            )
                    except:
                        for fastq in fastq_file:
                            if not fastq.endswith(".fastq.gz") and not fastq.endswith(".fq.gz"):
                                print_error(
                                    f"ERROR: FASTQ file does not have extension '.fastq.gz' or '.fq.gz': {fastq}",
                                    "Line",
                                    line,
                                )
                            if not exists(fastq):
                                print_error(
                                    f"ERROR: FASTQ file does not exist at this path: {fastq}",
                                    "Line",
                                    line,
                                )


main()
