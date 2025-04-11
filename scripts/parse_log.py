#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import re
from datetime import datetime
import pandas as pd

# patterns to match
time_pattern = r"\[\w{3} \w{3} \s?\d{1,2} \d{2}:\d{2}:\d{2} \d{4}\]"
rule_pattern = r"rule (\S+):"
jobid_pattern = r"jobid: (\d+)"

# setup the output table
df = pd.DataFrame(columns=["jobid", "rule", "start", "end"])
df.set_index("jobid", inplace=True)

# read the log file one line at a time
file_path = input("Path to log file: ")
with open(file_path, "r") as log_file:
    for line in log_file:
        
        # match the patterns
        time_match = re.search(time_pattern, line)
        rule_match = re.search(rule_pattern, line)
        jobid_match = re.search(jobid_pattern, line)

        # capture timestamp
        if time_match:
            timestamp = datetime.strptime(time_match.group().strip("[]"), "%a %b %d %H:%M:%S %Y")
        
        # capture rule
        if rule_match:
            rule = rule_match.group(1)
        
        # capture jobid (this will not work if a job has retries...)
        if jobid_match:
            jobid = int(jobid_match.group(1))
            if jobid not in df.index:
                df.loc[jobid] = [rule, timestamp, pd.NA]
            else:
                df.loc[jobid, "end"] = timestamp

# find the time difference
df.end = pd.to_datetime(df.end)
df["diff"] = (df.end - df.start).astype("timedelta64[s]")

# get the output
print(f"The longest time: {df['diff'].max()}")
out_path = input("Output directory: ")
df.to_csv(out_path + "/parse_log.csv")