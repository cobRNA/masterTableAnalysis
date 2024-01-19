#!/usr/bin/env python3
"""
 Read transcripts from masterTable and write tsv formatted output
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
    "endSupport",
    "spliced",
    "refCompare",
    "currentCompare",
    "sampleN",
    "Samples_metadata",
]

sys.stdout.write("\t".join(HEADER) + "\n")

# Create tsv line
for l in sys.stdin:
    line = l.replace('""', "NA").replace('"', "").replace(";", "").split()
    catalogues = line[13].split(",")
    for entry in catalogues:
        catalogue = entry.split(".")[0]
        output_line = [
            line[0],
            line[2],
            line[6],
            line[9],
            line[11],
            catalogue,
            line[15],
            line[17],
            line[19],
            line[21],
            line[23],
            line[25],
        ]
        sys.stdout.write("\t".join(output_line) + "\n")
