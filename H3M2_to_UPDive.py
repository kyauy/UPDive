#!/usr/bin/env python

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

if len(sys.argv) != 3 or "-h" in sys.argv or "--help" in sys.argv:
    print("python H3M2_to_UPDive.py <H3M2 output bed file> <UPDive MAD ROH reference file>")
    print("- First argument has to be the H3M2 bed file.")
    print("- Second argument has to be the UPDive MAD ROH size by chromosome reference.")
    print("Output directory will be the same directory as H3M2 bed file.")
    sys.exit()

print('H3M2_to_UPDive.py is launched.')
print('These files will be used for UPDive analysis: ' + str(sys.argv))

## Import files
input_file = sys.argv[1]
updive_ref = sys.argv[2]

## Function
def H3M2_to_UPDive(input_file, updive_ref):
    """
    python H3M2_to_UPDive.py <H3M2 output  bed file> <UPDive MAD ROH reference file>
    This is the main function that will process H3M2 output bed file into a UPDive report.
    - First argument has to be the H3M2 bed file.
    - Second argument has to be the UPDive MAD ROH size by chromosome reference.

    Output directory will be the same directory as H3M2 bed file.
    """
    # get sample name
    sample_abspath = os.path.abspath(input_file)
    sample_path = os.path.dirname(sample_abspath)
    sample_name_ext = os.path.basename(input_file)
    sample_name = os.path.splitext(sample_name_ext)[0]
    sample_path_table = sample_path + "/" + sample_name  + "_UPDive_table.txt"
    sample_path_report = sample_path + "/" + sample_name + "_UPDive_report.txt"

    # import file
    updive_ref = pd.read_csv(updive_ref, sep='\t', header=None)
    input_file = pd.read_csv(input_file, sep='\t', header=None)

    # process input file to get total ROH file by chromosome
    input_file['ROH_size'] = input_file.iloc[:,2] - input_file.iloc[:,1]
    input_file_ROHbychr = input_file.groupby([0], sort=False, as_index=False).sum()

    # exclude chrX
    input_file_ROHbychr = input_file_ROHbychr[0:22]

    # keep only chr column and total ROH size per chromosome
    input_file_ROHbychr = input_file_ROHbychr.loc[:,[0,'ROH_size']]

    # add UPDive MAD of ROH size by chromosome reference
    input_file_ROHbychr_UPDiveref = pd.concat([input_file_ROHbychr.reset_index(drop=True), updive_ref[1]], axis=1)

    # get MAD score as robustscale function (substract MAD and scale with MAD)
    input_file_ROHbychr_UPDiveref['MAD'] = (input_file_ROHbychr_UPDiveref.iloc[:,1] - input_file_ROHbychr_UPDiveref.iloc[:,2])  / input_file_ROHbychr_UPDiveref.iloc[:,2]

    # get clean output data format
    output_data = input_file_ROHbychr_UPDiveref[[0,"MAD"]]
    output_data.columns = [sample_name,'UPDive_MAD']

    print("Data processing is done for sample " + sample_name + ". UPDive table and reports will be written.")

    # produce UPDive table
    output_data.to_csv(sample_path_table, sep='\t', index=False)

    # produce UPDive report
    sys.stdout = open(sample_path_report,"w")
    ### header
    print("UPDive report")
    print("Sample: " + sample_name)
    ### recognize consanguinity
    if output_data.UPDive_MAD[output_data['UPDive_MAD'] >= 3].count() > 2 :
        print("Consanguinity suspected: the number of chromosome with MAD > 3 is " + str(output_data.UPDive_MAD[output_data['UPDive_MAD'] >= 3].count()) + ".")
        print("Interpretation of the following results may be wrong.")
    ### report significant iUPD events
    print("Suspected Uniparental Isodisomy (iUPD) events (MAD > 22):")
    if output_data.UPDive_MAD[output_data['UPDive_MAD'] >= 22].count() > 0:
        print(output_data[output_data['UPDive_MAD'] >= 22].to_string(index=False))
    else:
        print("No suspected iUPD.")
    sys.stdout.close()

# Run
H3M2_to_UPDive(input_file, updive_ref)
