# -*- coding: utf-8 -*-
import os
import pandas as pd

# ask for the location of the reads
#reads_loc = input("Enter file location: ")
reads_loc = "/fs/cbsuclarkfs1/storage/ede34/finch_lowpass_2022/XH-4057/20241025_LH00179_0133_A22HH77LT4"
out_loc = "/workdir/reg259/project"

# create a dictionary
data = {"name":[], "read1":[], "read2":[]}

# append to the dictionary
for folder in sorted(os.listdir(reads_loc)):
    if "Sample_XH" in folder:
        for file in sorted(os.listdir(f"{reads_loc}/{folder}")):
            if "R1" in file:
                name = f"03Dap{file.split('_')[0].split('XH-4057-')[-1]}"
                data["name"].append(name)
                data["read1"].append(f"{reads_loc}/{folder}/{file}")
            elif "R2" in file:
                data["read2"].append(f"{reads_loc}/{folder}/{file}")

# setup the dataframe
df = pd.DataFrame(data)
df.to_csv(out_loc + "/samples.csv", index=False)