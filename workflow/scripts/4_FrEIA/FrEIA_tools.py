#!/usr/bin/env python3
# Author: Norbert Moldovan
import sys
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.colors as mcolors
import csv

# This script bins the reads.
ChrLen = {"chr1": 248956422,
          "chr2": 242193529,
          "chr3": 198295559,
          "chr4": 190214555,
          "chr5": 181538259,
          "chr6": 170805979,
          "chr7": 159345973,
          "chr8": 145138636,
          "chr9": 138394717,
          "chr10": 133797422,
          "chr11": 135086622,
          "chr12": 133275309,
          "chr13": 114364328,
          "chr14": 107043718,
          "chr15": 101991189,
          "chr16": 90338345,
          "chr17": 83257441,
          "chr18": 80373285,
          "chr19": 58617616,
          "chr20": 64444167,
          "chr21": 46709983,
          "chr22": 50818468,
          "chrM": 16569
          }


def Binning(SampleDf, args):
    if args.chrom_sizes:
        # Read chromosome sizes.
        try:
            with open(args.chrom_sizes, mode="r") as infile:
                reader = csv.reader(infile,
                                    delimiter="\t")
                ChrLen = {rows[0]: int(rows[1]) for rows in reader}
        except:
            sys.exit("ERROR: No chromosomes length file found! Please check "
                     " if the file format is <chromosome name> <tab> <length>!")
    # Create bin edges for the input binsize.
    ChrLenSubset = {key: value for key, value in ChrLen.items()
                    if key in SampleDf["chr"].unique()}
    BinnedDf = SampleDf[SampleDf["chr"].isin(list(ChrLen.keys()))
                        ].reset_index(drop=True)

    if not args.binSize:
        args.binSize = 5000000
    for c in BinnedDf["chr"].unique():
        binEdge = range(1, ChrLenSubset.get(c)+args.binSize, args.binSize)
        # Binning and labeling the bins.
        binlab = list(c + "_" + str(l) for l in range(1, len(binEdge)))

        BinnedDf.loc[BinnedDf.chr == c, "bin"] = pd.cut(BinnedDf[BinnedDf.chr == c]["leftmost_nc_pos"],
                                                        binEdge,
                                                        labels=binlab)
    BinnedDf["bin"] = BinnedDf["bin"].astype("category")
    # Cast string data into categorical for
    # better memory usage and faster calculations.

    return BinnedDf


def DetectControl(sampTDF):
    # Get the control and the affected groups and create a group list.
    ControlGroup = [c for c in sampTDF["group"] if "*" in c]
    ControlSet = set(ControlGroup)

    if len(ControlSet) == 0:
        sys.exit("ERROR: Missing control group! Please type a '*' next "
                 " to every control group name!")

    elif len(ControlSet) > 1:
        sys.exit("".join(("ERROR: More than one control group:", ControlSet)))
    ControlGroup = ControlGroup[0][:-1]
    SampleGroups = sampTDF["group"].str.replace("*", "", regex=False).unique()
    # Set the control as last in the list.
    SampleGroups = np.setdiff1d(SampleGroups, ControlGroup)
    SampleGroups = np.append(SampleGroups, ControlGroup)

    return {"ControlGroup": ControlGroup,
            "SampleGroups": SampleGroups.tolist()}


def RegroupSamples(rgrPath, DataDf):
    try:
        rgrDf = pd.read_csv(rgrPath, delim_whitespace=True)
    except KeyError:
        rgrDf = pd.read_csv(rgrPath)

    rgrDf["group"] = rgrDf["group"].str.replace("*", "", regex=False)

    # Get sample names.Palette
    SampleL = rgrDf["sample_name"].unique()
    # Check if the current sample should be dropped.
    if DataDf["WhichSample"][0] not in SampleL:
        # Empty the dataframe but keep column headers.
        DataDf = pd.DataFrame(columns=DataDf.columns)
    else:
        # Set the new group for the sample.
        DataDf["WhichGroup"] = rgrDf[rgrDf["sample_name"
                                           ] == DataDf["WhichSample"
                                                       ][0]]["group"].values[0]
    return DataDf


def GeneratePalette(WhichGroup, ControlGroup):
    colors = sns.color_palette("Set1", len(WhichGroup)-1)
    # colors.append([mcolors.to_rgb("silver")])
    # Append silver color to colors for the ControlGroup.
    colors = np.append(colors, [mcolors.to_rgb("dimgray")], axis=0)
    Palette = dict(zip(WhichGroup, colors))
    return Palette


def CastDataTypes(data):
    # Cast object data into categorical,
    # and int32 for better memory usage and faster calculations.
    # print(data.info(memory_usage="deep"))
    data = data.reset_index(drop=True)
    for c in data.columns:
        if type(data[c][0]) == "object":
            data[c] = data[c].astype("category")
        elif type(data[c][0]) == str:
            data[c] = data[c].astype("category")
        elif type(data[c][0]) == "NA":
            data[c] = data[c].astype("category")
        elif type(data[c][0]) == int:
            data[c] = data[c].astype("int32")
        elif type(data[c][0]) == float:
            data[c] = data[c].astype("float32")
    # print(data.info(memory_usage="deep"))
    return data
