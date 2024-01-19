#!/usr/bin/env python3
"""
Read line from stdin, add to gene_id and transcript_id dot separated string 
containing chromosome.start.stop.strand to make ids unique. Written for bash pipeline.
"""

import sys


for l in sys.stdin:
    line = l.split()

    chromosome = line[0]
    start_pos = line[3]
    end_pos = line[4]
    strand = line[6]
    gene_id = line[9].replace('"', "").replace(";", "")
    transcript_id = line[11].replace('"', "").replace(";", "")

    # create output string
    output_line = (
        "\t".join(line[:9])  # with 'gene_id'
        + " "  # tab needed to separate . from transcipt_id
        + f'"{gene_id}.{chromosome}.{start_pos}.{end_pos}.{strand}";'
        + " "
        + line[10]
        + " "
        + f'"{gene_id}.{chromosome}.{start_pos}.{end_pos}.{strand}";'
        + "\n"
    )
    sys.stdout.write(output_line)  # return the line to stdout
