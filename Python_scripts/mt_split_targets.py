#!/usr/bin/env python3
"""
 Read transcripts from mT, split each target into new line and write tsv formatted output
 lines to stdout. Written for bash pipeline.
"""

import sys

# Write header to stdout
HEADER = [
    "chr",
    "type",
    "strand",
    "gene_id",
    "transcript_id",
    "target",
]

sys.stdout.write("\t".join(HEADER) + "\n")

# Create tsv line
for l in sys.stdin:
    line = l.replace('""', "NA").replace('"', "").replace(";", "").split()
    targets = line[13].split(",")
    for target in targets:
        output_line = [line[0], line[2], line[6], line[9], line[11], target]
        sys.stdout.write("\t".join(output_line) + "\n")