#!/usr/bin/env python3
"""
 Read target_ids from stdin and for each target_id find mT lines containing
 this target_id in target column.
 Reports transcript only once, even if it is linked to multiple targets
 but only within process. In the last step outputs from different processes are joined together,
 duplicated lines are removed and resulting transcripts are saved to separate file.
 Requires Python >= 3.12.1 for itertools batched.
"""

import sys
import os
import argparse
import warnings
from multiprocessing import Process
from multiprocessing import active_children

# Handle comandline arguments
parser = argparse.ArgumentParser(
    description="Find lncRNA related transcripts using taget_ids.",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
parser.add_argument(
    "mt_file_path",
    help="path to transcripts ICmT file",
    default=None,
)

parser.add_argument(
    "targets_path",
    help="path to file with target_ids",
    default=None,
)

parser.add_argument(
    "-p",
    "--processes",
    help="number of processes",
    type=int,
    default=1,
)

parser.add_argument(
    "-d",
    "--dir",
    help="output directory name",
    type=str,
    default="output",
)


args = parser.parse_args()
config = vars(args)


def split_data(data_set, n):
    from itertools import batched

    if n < 1:
        raise ValueError("n must be at least one")

    part_size = int(len(data_set) / n)
    splited_data = list(batched(data_set, part_size))

    if len(splited_data) > n:
        last_part = splited_data.pop()
        splited_data[-1] = splited_data[-1] + last_part  # its a tuple here

    return splited_data


def find_transcripts(targets, transcripts, output_filename):
    transcripts_set = set(transcripts)
    file = open(output_filename, "a", encoding="utf-8")
    for transcript in transcripts_set:
        for target_id in targets:
            if target_id in transcript:
                file.write(transcript + "\n")
                break  # leave loop to avoid duplicates
            else:
                pass  # continue searching

    file.close()


def join_output(list_of_filenames):
    joined_output = set()
    for file_name in list_of_filenames:
        with open(file_name, "r", encoding="utf-8") as file:
            joined_output.union(set(file.read().split(sep="\n")))

    # Write output to file
    with open(f"{config['dir']}/joined_output", "w", encoding="utf-8") as file:
        file.write("\n".join(joined_output))


if __name__ == "__main__":
    # Check whether the specified path exists or not
    if os.path.exists(config["dir"]):
        raise SystemExit(
            "Output directory already exists.\n\
                Submit new directory name using --dir \
                    or remove output directory from the previous run.\n\
                        Exitting program."
        )

    os.makedirs(config["dir"])

    # Open transcripts file containing only transcripts and load into memory
    with open(config["mt_file_path"], "r", encoding="utf-8") as file:
        transcripts = set(file.read().split(sep="\n"))

    # Open targets file containing only transcripts and load into memory
    with open(config["targets_path"], "r", encoding="utf-8") as file:
        targets = set(file.read().split())

    # Check Python version. At least 3.12 is required for itertools batched and paralellism
    if ((sys.version_info[0] < 3) or (sys.version_info[1] < 12)) and (
        config["processes"] > 1
    ):
        warnings.warn(
            "Warning........... Python >= 3.12.1 required for multiprocessing!\n\
                      The program continues to run, but using single process."
        )
        # Run only 1 process
        config["processes"] = 1
        # create list containing 1 set for compatibility
        transcripts_chunks = [transcripts]

    number_of_parts = config["processes"]
    transcripts_chunks = split_data(n=number_of_parts, data_set=transcripts)

    processes = []
    file_names = []
    for n in range(config["processes"]):
        processes.append(
            Process(
                target=find_transcripts,
                args=(targets, transcripts_chunks[n], f"{config['dir']}/output_{n}"),
            )
        )
        file_names.append(f"output_{n}")
    for process in processes:
        process.start()
    # get a list of all active child processes
    children = active_children()
    # report a count of active children
    print(f"Active Children Count: {len(children)}")
    # report each in turn
    for child in children:
        print(child)
    # wait for all processes to finish
    for t in processes:
        t.join()

    # Join output_files after all processes finish


print("Job completed!\n")
