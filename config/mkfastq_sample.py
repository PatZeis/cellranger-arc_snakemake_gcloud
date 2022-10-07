#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 15:50:02 2022

@author: patrice.zeis
"""

import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--sample_ids", help="name of samples", type=str)
parser.add_argument("-l", "--lanes", help="give lanes", type=str)
parser.add_argument("-i", "--index", help="give index", type=str)
parser.add_argument("-f", "--flow_cells", help="give flowcells numbers(s)", type=str)
parser.add_argument("-n", "--flow_cell_names", help="give flowcell RNA name(s)", type=str)
parser.add_argument("-m", "--multiomic_atac_index",nargs='?' , help="give_atac_index", type=str)
parser.add_argument("-a", "--flow_cell_names_atac",nargs='?', help="give flowcell ATAC name(s)", type=str)

args = parser.parse_args()

flow_cell_names = args.flow_cell_names
flow_cells = args.flow_cells
lanes = args.lanes
samples = args.sample_ids
index = args.index
if args.multiomic_atac_index is not None:
    atac_index = args.multiomic_atac_index
    atac_index = atac_index.split(".") 
if args.flow_cell_names_atac is not None:
    flow_cell_names_atac = args.flow_cell_names_atac
    flow_cell_names_atac = flow_cell_names_atac.split(".")

samples = samples.split(".")

flow_cells = flow_cells.split(".")

lanes = lanes.split(".")

flow_cell_names = flow_cell_names.split(".")

index = index.split(".")

#rna1 = "R" in layer
#rna2 = "r" in layer
#rna = rna1 | rna2

#atac1 = "A" in layer[0]
#atac2 = "a" in layer[0]
#atac = atac1 | atac2

#if rna | atac:
#
#    if rna:
#        layer = "RNA"
#    else:
#        layer = "ATAC"
#
d = {}

for s, i in zip(samples, index):
    index_split = i.split(",")
    if len(index_split) == 2:
        index2=index_split[1]
        i = index_split[0]
    else:
        index2="."
    d[s] = [s,i, index2]

if len(flow_cell_names) == 1:
    flow_cell_names = "".join(flow_cell_names)
    for s in samples:
        d[s].append(flow_cell_names)
else:
    for s, m in zip(samples, flow_cell_names):
        d[s].append(m)

if len(flow_cells) == 1:
    flow_cells = "".join(flow_cells)
    for s in samples:
        d[s].append(flow_cells)
else:
    for s, f in zip(samples, flow_cells):
        d[s].append(f)


if len(lanes) == 1:
    lanes = "".join(lanes)
    lanes = ",".join(list(lanes))
    for s in samples:
        d[s].append(lanes)
else:
    for s, l in zip(samples, lanes):
        l =  ",".join(list(l))
        d[s].append(l)

if args.flow_cell_names_atac is not None:
    if len(flow_cell_names_atac) ==1:
        flow_cell_names_atac = "".join(flow_cell_names_atac)
        for s in samples:
            d[s].append(flow_cell_names_atac)
    else:
        for s, a in zip(samples, flow_cell_names_atac):
            d[s].append(a)

    if args.multiomic_atac_index is not None:
        for s, m in zip(samples, atac_index):
            d[s].append(m)
    else:
        raise ValueError("Please also give atac indixes")

if args.flow_cell_names_atac is not None:
    data = pd.DataFrame.from_dict(d, orient='index', columns=["id","index","index2","flowcell_name","flowcell_number","lane_numbers","flowcell_name_atac","atac_index"])
else:
    data = pd.DataFrame.from_dict(d, orient='index', columns=["id","index","index2","flowcell_name","flowcell_number","lane_numbers"])

data.to_csv("mkfastq_samples.tsv",sep='\t', index=False)
    
#else:
#    raise ValueError("If Atac layer give either 'ATAC', 'A', ','[Aa]tac' or 'a' and if RNA give either 'RNA', 'R', '[Rrna]' and 'r'")
