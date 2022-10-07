#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 08:12:00 2022

@author: patrice.zeis
"""

import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--RNA_path", help="give path to ATAC Fastqs", type=str)
parser.add_argument("--ATAC_path", help="give path to RNA Fastqs", type=str)
parser.add_argument("--wildcard", help="give sample", type=str)
parser.add_argument("--output_file", help="give output file name", type=str)
parser.add_argument("--samples", help="give input samples file", type=str)


args = parser.parse_args()
sample_file = args.samples
rna_path = args.RNA_path
atac_path = args.ATAC_path
wildcard = args.wildcard
output_file = args.output_file

samples_mk = pd.read_csv(sample_file, sep="\t", dtype=object)

def get_entry_mk(sample_id, column_name):
    entry = samples_mk.loc[samples_mk["id"] == sample_id, column_name].values
    if len(entry) > 1:
        raise ValueError("Unexpected number of items. Sample IDs are not unique")
    elif len(entry) < 1:
        raise ValueError( f"Unexpected number of items. Sample ID {sample_id} not found for columns {column_name}" )
    return entry[0]

flow_cell = get_entry_mk(wildcard, "flowcell_number")

sample_nam = wildcard+"_"+flow_cell

d = {}

d["rna"] = [rna_path, sample_nam, "Gene Expression"]
d["atac"] = [atac_path, sample_nam, "Chromatin Accessibility"]

data = pd.DataFrame.from_dict(d, orient='index', columns=["fastqs","sample","library_type"])
 
data.to_csv(output_file,sep=',', index=False, header=True)
