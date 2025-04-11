#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import pandas as pd

# change to the correct directory
workdir = "/workdir/ede34/projects/galap_finches/reference"
os.chdir(workdir)

# identify files
old_file = [file for file in os.listdir(workdir) if ".fna" in file][0]
new_file = "small_tree_finch_genome.fna"

# find all the scaffold names
info_file = [file for file in os.listdir(workdir) if ".tsv" in file][0]
info = pd.read_csv(info_file, sep="\t")

# open the original fasta file and create a new fasta file
with open(old_file, "r") as in_ref, open(new_file, "w") as out_ref:
        for line in in_ref:
            # if the line is a header, change it
            if line.startswith(">"):
                mask = info["GenBank seq accession"].apply(lambda x: x in line)
                new_line = ">" + info[mask]["Sequence name"].values[0] + "\n"
                out_ref.write(new_line)
            # otherwise keep it the same
            else:
                out_ref.write(line)
