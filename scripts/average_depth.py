# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np

df = pd.DataFrame(columns=["name", "chr", "count", "mean", "std", "min", "25%", "50%", "75%", "max"])

for file in snakemake.input:
    # extract values from each file
    data = pd.read_csv(file, sep="\t", header=None)
    data.rename(columns={0:"chromosome", 1:"position", 2:"depth"}, inplace=True)
    
    # append to the dataframe
    name = file.split("/")[-2]
    chr_num = file.split("/")[-1].split(".")[0]
    df.loc[df.shape[0]] = [name] + [chr_num] + data["depth"].describe().tolist()

# find the weighted averages
weighted_avg = lambda x: np.average(x["mean"], weights=x["count"])
df_sample = df.groupby("name").apply(weighted_avg, include_groups=False).reset_index(name="wt_avg")
df_chr = df.groupby("chr").apply(weighted_avg, include_groups=False).reset_index(name="wt_avg")

# create the output files
df.to_csv(snakemake.output[0], index=False)
df_sample.to_csv(snakemake.output[1], index=False)
df_chr.to_csv(snakemake.output[2], index=False)