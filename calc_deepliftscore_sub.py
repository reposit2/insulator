#!/home/user/.pyenv/shims/python3

# This code was created based on codes in 'Tasaki, S., Gaiteri, C., Mostafavi, S. & Wang, Y. Deep learning decodes the principles of differential gene expression. Nature Machine Intelligence (2020)'

import sys
args = sys.argv
print(args[0])
print(args[1])
print(args[2])
print(args[3])
print(args[4])
print(args[5])

deg_data_file = "./data/fourcelltypeswithhff_fpkm.txt.gz"
promoter_data_loc = "./data/Promoter_features/"
promoter_annotation_data = ["GTRDC"]
enhancer_data_loc2c = "./data/Promoter_features/"
enhancer_annotation_data2c = ["GTRDE4PROM2006"]
enhancer_data_loc2 = args[1]
enhancer_annotation_data2 = [args[2]]
test_genes2 = args[3]
outloc=args[4]
outfile=args[5]

import os
sys.path.append('./functions/')
import data
import model_utils
import pandas as pd
import numpy as np
import json
import subprocess
import time
import tensorflow as tf
import gzip

# Prepare data
shuffle="None"
Y_test2, X_promoter_test2, X_enhancer_test2 = data.prep_ml_data_split3(
    deg_data_file=deg_data_file,
    mRNA_data_loc=promoter_data_loc,
    mRNA_annotation_data=promoter_annotation_data,
    promoter_data_loc=enhancer_data_loc2,
    promoter_annotation_data=enhancer_annotation_data2,
    test_genes=test_genes2,
    outloc=outloc,
    shuffle=shuffle)

Y_test3, X_promoter_test3, X_enhancer_test3 = data.prep_ml_data_split3(
    deg_data_file=deg_data_file,
    mRNA_data_loc=promoter_data_loc,
    mRNA_annotation_data=promoter_annotation_data,
    promoter_data_loc=enhancer_data_loc2c,
    promoter_annotation_data=enhancer_annotation_data2c,
    test_genes=test_genes2,
    outloc=outloc,
    shuffle=shuffle)

# Compute Spearman's correlation between actual and predicted expression for each sample
print('Compute Spearmans correlation between actual and predicted expression for each sample')
# Read the best model
rsc = f"cat {outloc}summary/best_model.txt"
print(rsc)
subprocess.call(rsc, shell=True)
while(not os.path.isfile(f"{outloc}summary/best_model.txt")):
    time.sleep(1)
with open(outloc+"summary/best_model.txt") as f:
    best_model=f.readline().rstrip()
    print(best_model)

# Prediction for test samples with the best model
model_utils.test_prediction2(outloc,
                            best_model,
                            X_enhancer_test2,
                            Y_test2)
# Descriptions of outputs
#.{outloc}/test_data/prediction.txt.gz: Predicted gene expression data
#.{outloc}/test_data/actual.txt.gz: Actual gene expression data
#.{outloc}/test_data/geneid.txt.gz: Genes in testing data.

rc = f"{outloc}/test_data/geneid.txt.gz"
with gzip.open(rc) as f:
    geneidfile = f.read()
geneidfile = geneidfile.decode()
geneidlist = geneidfile.split('\n')
geneidlist2 = [a for a in geneidlist if a != '']
if len(geneidlist2) == 1:
    quit()

# Compute prediction accuracy for each sample
rsc = f"Rscript --vanilla --slave functions/calc_performance.R {outloc} &> /dev/null"
print(rsc)
try:
    subprocess.call(rsc, shell=True)
except:
    quit()
while(not os.path.isfile(f"{outloc}test_data/cor_tbl.txt")):
    time.sleep(1)
# Descriptions of outputs
#.{outloc}/test_data/cor_tbl.txt: Correlation between acutual and predicted gene expression for each sample
rc = f"{outloc}test_data/cor_tbl.txt"
print(rc)
pd.read_csv(rc,sep="\t")

# Compute average DeepLIFT score for each regulator
print('Compute average DeepLIFT score for each regulator')
# Read the best model
with open(outloc+"summary/best_model.txt") as f:
    best_model=f.readline().rstrip()

# Estimate variable imporance using test samples
model_utils.compute_DeepLIFT4(outloc,
                             best_model,
                             X_enhancer_test2,
                             Y_test2,
                             X_enhancer_test3,
                             Y_test3)

# Descriptions of outputs
#.{outloc}/DeepLIFT/DNA_{sample index}.txt.gz: DeepLIFT scores of promoter and enhancer regulators for each sample. The sample index corresponds to the column index of gene expression data, which starts from 0. DeepLIFT scores of regulators were sparated with commas. Its order is identical to the one appeared in {outloc}/feature_norm_stats.txt

# Concatenate and summarize DeepLIFT scores 
rsc = f"Rscript --vanilla --slave functions/summarize_DeepLIFT.R {outloc} &> /dev/null"
print(rsc)
subprocess.call(rsc, shell=True)
while(not os.path.isfile(f"{outloc}DeepLIFT/enhancer_importance_mean.txt")):
    time.sleep(1)
# Descriptions of outputs
#.{outloc}/DeepLIFT/enhancer_importance_mean.txt: Average DeepLIFT score for each promoter and enhancer regulator in each tissue.
rc = f"{outloc}DeepLIFT/enhancer_importance_mean.txt"
print(rc)
pd.read_csv(rc,sep="\t")
with open(outfile, 'a') as f:
    print(rc, file=f)

