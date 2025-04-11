#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import pandas as pd

# location of the reads and the records sheet
reads_loc = "/fs/cbsuclarkfs1/storage/ede34/finch_lowpass_2022/XH-4057/20241025_LH00179_0133_A22HH77LT4"
out_loc = "/home2/reg259/finch_genome_imputation/example_project"
records = pd.read_csv("/home2/reg259/finch_genome_imputation/example_project/Daphne_finches_band_records_2022_CLEANED.csv")

# create a dictionary
data = {"name":[], "read1":[], "read2":[]}

# look in every folder for possible samples
for folder in sorted(os.listdir(reads_loc)):
    if "Sample_XH" in folder:
        
        # check if the sample is in the records
        num = folder.split("-")[2]
        if records[records["NUM.07"]==num].size:
            
            # choose the first meaningful unique value if there are multiple
            name = records[records["NUM.07"]==num]["meaningful.unique"].values[0]
            # need to deal with the one that was sequenced twice
            if len(folder.split("-")) > 3:
                name = name + "-" + folder.split("-")[3]
            data["name"].append(name)
            
            # append the reads to the dictionary
            for file in sorted(os.listdir(f"{reads_loc}/{folder}")):
                if "R1" in file:
                    data["read1"].append(f"{reads_loc}/{folder}/{file}")
                elif "R2" in file:
                    data["read2"].append(f"{reads_loc}/{folder}/{file}")

# setup the dataframe
df = pd.DataFrame(data)
df.to_csv(out_loc + "/samples2.csv", index=False)