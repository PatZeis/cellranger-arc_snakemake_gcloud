#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  6 08:05:44 2022

@author: patrice.zeis
"""

import argparse
import pandas as pd
import yaml

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_sample_file", help="give mkfastq sample file", type=str)
parser.add_argument("-b", "--bucket", help="give bucket name", type=str)

args = parser.parse_args()

bucket = args.bucket
sample_file = args.input_sample_file

sample_file = pd.read_csv(sample_file, sep='\t')

flowcell_name = sample_file.loc[:,"flowcell_name"][0] ### assuming that all samples are run on the same flowcell
flowcell_name_atac = sample_file.loc[:,"flowcell_name_atac"][0]
lanes = list(sample_file["lane_numbers"].to_numpy())[0]
samples = list(sample_file["id"].to_numpy())
flowcell = str(list(sample_file["flowcell_number"])[0])

d = {}
d["lanes"] = lanes.split(",")
d["fc"] = flowcell.split(";")
d["gs_bucket"] = bucket
d["references"] = 'config/references.tsv'
d["filter_gtf"] = False

d["bcl_file"] = 'data/mkfastq/bcl_folder'

d["reference_files"] = ["Genome", "SAindex", "chrNameLength.txt", "exonInfo.tab", "sjdbInfo.txt", "transcriptInfo.tab", "Log.out", "chrLength.txt", "chrStart.txt", 
                        "geneInfo.tab", "sjdbList.fromGTF.out.tab", "SA", "chrName.txt", "exonGeTrInfo.tab", "genomeParameters.txt", "sjdbList.out.tab"]

d["star_files"] = ['Aligned.sortedByCoord.out.bam','Log.final.out', 'Log.progress.out', 'SJ.out.tab', 'Solo.out/Barcodes.stats', 'Solo.out/Gene/Features.stats', 'Solo.out/Gene/Summary.csv','Solo.out/Gene/UMIperCellSorted.txt','Solo.out/Gene/raw/barcodes.tsv', 'Solo.out/Gene/raw/features.tsv','Solo.out/Gene/raw/matrix.mtx', 'Solo.out/Gene/filtered/barcodes.tsv','Solo.out/Gene/filtered/features.tsv','Solo.out/Gene/filtered/matrix.mtx','Solo.out/GeneFull/Features.stats', 'Solo.out/GeneFull/Summary.csv', 'Solo.out/GeneFull/UMIperCellSorted.txt','Solo.out/GeneFull/raw/barcodes.tsv','Solo.out/GeneFull/raw/features.tsv',  'Solo.out/GeneFull/raw/matrix.mtx', 'Solo.out/GeneFull/filtered/barcodes.tsv', 'Solo.out/GeneFull/filtered/features.tsv','Solo.out/GeneFull/filtered/matrix.mtx',  'Solo.out/Velocyto/Features.stats',  'Solo.out/Velocyto/Summary.csv','Solo.out/Velocyto/raw/barcodes.tsv', 'Solo.out/Velocyto/raw/features.tsv', 'Solo.out/Velocyto/raw/matrix.mtx'    ]

d["read"] = ['R1','R2']

d["RNA_mkfastq_reads"] = ['R1','R2','I1','I2']

all_html = 'outs/fastq_path/Reports/html/'+flowcell_name+'/all/all/all/lane.html'
all_laneBarcode_html = 'outs/fastq_path/Reports/html/'+flowcell_name+'/all/all/all/laneBarcode.html'
all_default_html = 'outs/fastq_path/Reports/html/'+flowcell_name+'/default/all/all/lane.html'
all_default_laneBarcode_html = 'outs/fastq_path/Reports/html/'+flowcell_name+'/default/all/all/laneBarcode.html'
all_default_undetermined_html = 'outs/fastq_path/Reports/html/'+flowcell_name+'/default/Undetermined/all/lane.html'
all_default_undetermined_laneBarcode_html = 'outs/fastq_path/Reports/html/'+flowcell_name+'/default/Undetermined/all/laneBarcode.html'
default_undertermined_unknown_html = 'outs/fastq_path/Reports/html/'+flowcell_name+'/default/Undetermined/unknown/lane.html'
default_undertermined_unknown_laneBarcode_html = 'outs/fastq_path/Reports/html/'+flowcell_name+'/default/Undetermined/unknown/laneBarcode.html'
flowcell_name_all_all_html = 'outs/fastq_path/Reports/html/'+flowcell_name+'/'+flowcell_name+'/all/all/lane.html'
flowcell_name_all_all_laneBarcode_html = 'outs/fastq_path/Reports/html/'+flowcell_name+'/'+flowcell_name+'/all/all/laneBarcode.html'

mkfastq_files = ['outs/input_samplesheet.csv','outs/interop_path/IndexMetricsOut.bin', 'outs/fastq_path/Reports/html/Report.css', 'outs/fastq_path/Reports/html/index.html', 'outs/fastq_path/Reports/html/tree.html',
                 all_html, all_laneBarcode_html, all_default_html, all_default_laneBarcode_html, all_default_undetermined_html , all_default_undetermined_laneBarcode_html, default_undertermined_unknown_html,default_undertermined_unknown_laneBarcode_html,
                 flowcell_name_all_all_html, flowcell_name_all_all_laneBarcode_html  ]

d["mkfastq_files"] = mkfastq_files


d["fastq_path"] = "outs/fastq_path/"+flowcell_name
d["fastq_path_atac"] = "outs/fastq_path/"+flowcell_name_atac
d["fastq_undetermined"] = "outs/fastq_path"
d["flow_cell"] = flowcell_name
d["flow_cell_atac"] = flowcell_name_atac

with open('config.yaml', 'w') as file:
    documents = yaml.dump(d, file, default_flow_style=False)

        
        




