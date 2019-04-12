########
## H3M2 to UPDive report, via MAD normalization
## Author : Kevin Yauy
## Input : bed file from H3M2
## Date : 2019-04-12
########

## Library
import os
import sys
import numpy as np
import pandas as pd

## OS interface
print('H3M2 to UPDive report is launched.')
print('This file will be use for processing: ', str(sys.argv))

## Import files
input_file = sys.argv[0]
updive_ref = sys.argv[1]

## Function
def H3M2_to_UPDive(input_file, updive_ref):
    # get sample name
    sample_path = os.path.basename(input_file)
    sample_name = os.path.splitext(sample_path)[0]
    sample_path_table = sample_path  + "_UPDive_table.txt"
    sample_path_report = sample_path + "_UPDive_report.txt"
    print(sample_path_report)

    # import file
    updive_ref = pd.read_csv(updive_ref, sep='\t', header=None)
    print(updive_ref.head())
    input_file = pd.read_csv(input_file, sep='\t', header=None)
    print(input_file.head())

    # process input file to get total ROH file by chromosome
    input_file['ROH_size'] = input_file[2] - input_file[1]
    input_file_ROHbychr = input_file.groupby([0], sort=False, as_index=False).sum()
    print(input_file_ROHbychr.head())

    # exclude chrX
    input_file_ROHbychr = input_file_ROHbychr[0:22]
    print(input_file_ROHbychr.head())

    # keep only chr column and total ROH size per chromosome
    input_file_ROHbychr = input_file_ROHbychr.loc[:,[0,'ROH_size']]
    print(input_file_ROHbychr.head())

    # add UPDive MAD of ROH size by chromosome reference
    input_file_ROHbychr_UPDiveref = pd.concat([input_file_ROHbychr.reset_index(drop=True), updive_ref[1]], axis=1)
    print(input_file_ROHbychr_UPDiveref.head())

    # get MAD score as robustscale function (substract MAD and scale with MAD)
    input_file_ROHbychr_UPDiveref['MAD'] = (input_file_ROHbychr_UPDiveref.iloc[:,1] - input_file_ROHbychr_UPDiveref.iloc[:,2])  / input_file_ROHbychr_UPDiveref.iloc[:,2]

    # get clean output data format
    output_data = input_file_ROHbychr_UPDiveref[[0,"MAD"]]
    output_data.columns = [sample_name,'UPDive_MAD']

    output_data.to_csv(sample_path_table, sep='\t')

    sys.stdout = open(sample_path_report,"w")

    # produce UPDive report
    ### header
    print("UPDive report")
    print("Sample: " + sample_name)
    ### recognize consanguinity
    if output_data.UPDive_MAD[output_data['UPDive_MAD'] >= 3].count() > 2 :
        print("Consanguinity suspected: the number of chromosome with MAD > 3 is " + str(output_data.UPDive_MAD[output_data['UPDive_MAD'] >= 3].count()) + ".")
        print("Interpretation of the following results may be wrong.")

    print("Suspected UPD events:")
    print(output_data.UPDive_MAD[output_data['UPDive_MAD'] >= 22])
    sys.stdout.close()

# Run
H3M2_to_UPDive(input_file, updive_ref)

