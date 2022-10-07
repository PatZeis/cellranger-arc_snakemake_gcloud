#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 10:03:27 2022

@author: patrice.zeis
"""

import pandas as pd
import numpy as np
import argparse



parser = argparse.ArgumentParser()
parser.add_argument("--samples", help="give mkfastq_samples file", type=str)
parser.add_argument("--wildcard", help="give sample", type=str)
parser.add_argument("--output_file", help="give output file name", type=str)
parser.add_argument("--atac", nargs='?' , help="create atac simple csv file from mkfastq_samples file", type=str)


args = parser.parse_args()

sample_file = args.samples
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
lanes = get_entry_mk(wildcard, "lane_numbers")
sample = get_entry_mk(wildcard,"id")
index1 = get_entry_mk(wildcard,"index") 
index2 = get_entry_mk(wildcard,"index2")

if args.atac is None:

    if index2 == ".":
        index = index1
    else:
        index = [index1, index2]
        index = ",".join(index)

    lanes = lanes.split(",")
    index = index.split(",")

    sample_nam = sample+"_"+flow_cell


    d = {}
    if len(index) == 1:
        index = "".join(index)
        for l in lanes:

            lanes_id =  [l,sample_nam, index]
            d[l] = lanes_id
    
    elif len(index) == 2:
        for l in lanes:
            index1 = index[0]
            index2 = index[1]
   
            lanes_id =  [l,sample_nam,index1,index2]
            d[l] = lanes_id
        
    else:
        raise ValueError("Give either 1 or 2 index, something went wrong for the mkfastq_sample file")
    
    if len(d[l]) == 4:
        data = pd.DataFrame.from_dict(d, orient='index', columns=["Lane","Sample","Index","Index2"])

    else:
        data = pd.DataFrame.from_dict(d, orient='index', columns=["Lane","Sample","Index"])
else:
    index = get_entry_mk(wildcard,"atac_index") 
    sample_nam = sample+"_"+flow_cell
    lanes = lanes.split(",")
    d = {}
    for l in lanes:
        lanes_id =  [l,sample_nam, index]
        d[l] = lanes_id
    
    data = pd.DataFrame.from_dict(d, orient='index', columns=["Lane","Sample","Index"])


data.to_csv(output_file,sep=',', index=False, header=True)
