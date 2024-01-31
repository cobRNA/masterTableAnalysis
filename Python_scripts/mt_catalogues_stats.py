#!/usr/bin/env python3
"""
 Read transcripts from primary_target file and write tsv formatted output
 lines to stdout. Written for bash pipeline.
"""

import sys

# Write header to stdout
HEADER = [
    "chr",
    "CLS",
    "type",
    "start",
    "end",
    "strand",
    "target_id",
    "target_set",
    "target_category",
]

sys.stdout.write("\t".join(HEADER) + "\n")

# Create tsv line
for l in sys.stdin:
    line = l.replace('""', "NA").replace('"', "").replace(";", "").split()
    output_line = [*line[0:5], line[6], line[9], line[11], line[13]]
    sys.stdout.write("\t".join(output_line) + "\n")
