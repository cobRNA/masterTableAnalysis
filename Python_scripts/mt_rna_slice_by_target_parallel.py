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


def concat_output(list_of_filenames):
    concat_output = set()
    concat_file_name = "concatenated_output.txt"
    for file_name in list_of_filenames:
        with open(file_name, "r", encoding="utf-8") as file:
            concat_output.union(set(file.read().split(sep="\n")))

    # Write output to file
    with open(
        os.path.join(config["dir"], concat_file_name), "w", encoding="utf-8"
    ) as file:
        file.write("\n".join(concat_output))

    # Repot status
    print("Completed.")
    print(
        f"Concatenated output file location: {os.path.join(config['dir'], concat_file_name)}"
    )


if __name__ == "__main__":
    # Check whether the specified path exists or not
    if os.path.exists(config["dir"]):
        raise SystemExit(
            """Output directory already exists.\nSubmit new directory name using --dir or remove output directory from the previous run.\nExitting program."""
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
        print("Python >= 3.12.1 required for multiprocessing!")
        print("The program continues to run, but using single process.")

        # Run only 1 process
        config["processes"] = 1
        # create list containing 1 set for compatibility
        transcripts_chunks = [transcripts]

    # Fix negative process number, do not split the data if only 1 process is going to run
    elif config["processes"] <= 1:
        # Run only 1 process
        config["processes"] = 1
        # create list containing 1 set for compatibility
        transcripts_chunks = [transcripts]

    # If requirements are met and more processes required, split data using itertools
    else:
        number_of_parts = config["processes"]
        transcripts_chunks = split_data(n=number_of_parts, data_set=transcripts)

    processes = []
    file_names = []
    for n in range(config["processes"]):
        processes.append(
            Process(
                target=find_transcripts,
                args=(
                    targets,
                    transcripts_chunks[n],
                    os.path.join(config["dir"], f"output_{n}"),
                ),
            )
        )
        file_names.append(os.path.join(config["dir"], f"output_{n}"))
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
    for p in processes:
        p.join()

    # Report status
    print("Search completed!")

    # Concatenate output_files after all processes finish
    if len(file_names) > 1:
        print("Concatenating output files...")
        concat_output(file_names)

    # Report job status
    print("~~~ Job completed! ~~~")
