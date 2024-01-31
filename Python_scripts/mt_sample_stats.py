#!/usr/bin/env python3

"""
 1. Read transcript from stdin;
 2. Find transcript_id and sample metadata code/s;
 3. Split samples into different lines;
 4. Substitude sample codes with sample metadata from Hv3_sampleCoding file;
 5. Return transcript_id \t full sample metadata line into stdout.
"""


import sys

# Create dictionary to store metadata
meta_data = {}
# Read metadata from file
PATH = "/home/tomasz/Desktop/mT/Hv3_sampleCoding"
with open(PATH, "r", encoding="utf-8") as file:
    lines = file.readlines()
    for l in lines:
        line = l.split()
        meta_data[line[0]] = "\t".join(line[1:])


for l in sys.stdin:
    line = l.split()
    codes = line[25].strip('";').split(",")
    transcript_id = line[11].strip('";')
    for code in codes:
        # create output string
        output_line = transcript_id + "\t" + meta_data[code] + "\n"
        sys.stdout.write(output_line)  # return the line to stdout