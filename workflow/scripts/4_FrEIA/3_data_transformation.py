#!/usr/bin/env python3
# Norbert Moldovan

import os
import sys
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from FrEIA_tools import Binning, CastDataTypes


# Parsing the arguments + some nice help text.
def ParsingArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input",
                        dest="inPath",
                        type=str,
                        required=True,
                        help="The "
                        "path to the sample file produced in extraction step.")
    parser.add_argument("-o", "--output",
                        dest="output",
                        type=str,
                        required=True,
                        help="The path to the outputput folder.")
    parser.add_argument("-st", "--sampTable",
                        dest="SampleTable",
                        type=str,
                        required=True,
                        help="The path to the original sample sheet.")
    parser.add_argument("-fsmin", "--fragmSizeMin",
                        dest="fraLenMin",
                        type=int,
                        required=False,
                        help="The minimal fragment size.")
    parser.add_argument("-fsmax", "--fragmSizeMax",
                        dest="fraLenMax",
                        type=int,
                        required=False,
                        help="The maximal fragment size.")
    parser.add_argument("-subs", "--subsample",
                        dest="Subsample",
                        type=float,
                        required=False,
                        help="Subsample the data randomly. Example: "
                        "'0.5' - returns a randomly sampled dataframe "
                        "50%% of the original size.")
    parser.add_argument("--bootstrap_sample",
                        dest="SubSampRepeat",
                        type=int,
                        required=False,
                        help="Bootstrap sample the reads and creat this many "
                        "outputs with randomly sampled reads from the input. "
                        "Same reads can be included multiple times.")
    parser.add_argument("-bs", "--binSize",
                        dest="binSize",
                        type=int,
                        required=False,
                        help="The size of the genomic bins for "
                        "subchromosomal analysis. The default is 5,000,000 bp")
    parser.add_argument("-b", "--both",
                        dest="doBoth",
                        action="store_true",
                        help="Calculate statistics for both read ends.")
    parser.add_argument("-schoff", "--subchroff",
                        dest="subchroff",
                        action="store_false",
                        default="True",
                        help="Turn off subchromosomal analysis. "
                        "Subchromosomal analysis can use a lot of memory"
                        "[>16Gb] for large genomes. Default is On.")
    return parser.parse_args()


Cols = ["qname",
        "chr",
        "read_len",
        "strand",
        "leftmost_nc_pos",
        "rightmost_nc_pos",
        "leftmost_nc_seq",
        "leftmost_tnc_seq",
        "rightmost_tnc_seq",
        "rightmost_nc_seq",
        "leftmost_read_seq",
        "rightmost_read_seq",
        "gccont"
        ]

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
          "chrM": 16569}

"""=== Dataset transformations ==="""
# Calculate the proportion of features by length, genomic bin,
# chromosome and sample.
# Proportions by length per chromosome are also calculated!


def BaseCompAbCalc(SampleDf, group, whichEnd, lvl):
    BCDf = pd.DataFrame(columns=["Unit", "Base", "Length", "CountSamp",
                                 "CountLen", "RelAbSamp", "WhichEnd",
                                 "WhichGroup"])
    OutDf = pd.DataFrame()

    for e in whichEnd:
        # print(SampleDf[e][0], len(SampleDf))
        if lvl == "Genome_Lvl":
            hBCDf = SampleDf.groupby(["read_len", e], observed=True
                                     ).strand.count().reset_index()
            BCDf.CountLen = hBCDf.strand
            BCDf.CountSamp = hBCDf.groupby([e], observed=True
                                           ).strand.transform("sum")
            BCDf.RelAbSamp = (hBCDf.groupby([e], observed=True
                                            ).strand.transform("sum")
                                                / len(SampleDf.chr))
            BCDf.Unit = "GRCh38"

        elif lvl == "Chromosome_Lvl":
            hBCDf = SampleDf.groupby(["chr", "read_len", e], observed=True
                                     ).strand.count().reset_index()
            BCDf.CountLen = hBCDf.strand
            BCDf.CountSamp = hBCDf.groupby(["chr", e], observed=True
                                           ).strand.transform("sum")
            BCDf.RelAbSamp = (hBCDf.groupby(["chr", e], observed=True
                                            ).strand.transform("sum")
                                    / hBCDf.groupby(["chr"], observed=True
                                                    ).strand.transform("sum"))
            BCDf.Unit = hBCDf.chr

        elif lvl == "SubChr_Lvl":
            hBCDf = SampleDf.groupby(["bin", "read_len", e], observed=True
                                     ).strand.count().reset_index()
            BCDf.CountLen = hBCDf.strand
            BCDf.CountSamp = hBCDf.groupby(["bin", e], observed=True
                                           ).strand.transform("sum")
            BCDf.RelAbSamp = (hBCDf.groupby(["bin", e], observed=True
                                            ).strand.transform("sum")
                                    / hBCDf.groupby(["bin"], observed=True
                                                    ).strand.transform("sum"))
            BCDf.Unit = hBCDf.bin

        BCDf.Base = hBCDf[e]
        BCDf.Length = hBCDf.read_len
        BCDf.WhichEnd = e.split("_")[0]
        BCDf.WhichGroup = group

        # Cast string data into categorical for better memory usage and
        # faster calculations.
        for c in ["Unit", "Base", "WhichEnd", "WhichGroup"]:
            BCDf[c] = BCDf[c].astype("category")

        OutDf = OutDf.append([BCDf]).reset_index(drop=True)
        # print(OutDf, group, lvl)
        # print(OutDf.info(memory_usage="deep"), "\n")
    return OutDf


# Mononucleotide data is processed here.
def MonoNcAnalysis(SampleDf, args, filename, group, lvl):
    # Set ends needed to be analyzed.
    if args.doBoth:
        whichEnd = ["leftmost_nc_seq", "rightmost_nc_seq"]
    else:
        whichEnd = ["leftmost_nc_seq"]

    BCDf = BaseCompAbCalc(SampleDf, group, whichEnd, lvl)

    # Save the results.
    BCDf.to_parquet("".join((args.output, "4_FrEIA/3_Abundances/",
                             group, "/", lvl, "/M__", filename, ".pq")),
                    compression="gzip",
                    index=False)


# Trinucleotide data is processed here.
def TriNcAnalysis(SampleDf, args, filename, group, lvl):

    # Set ends needed to be analyzed.
    if args.doBoth:
        whichEnd = ["leftmost_tnc_seq", "rightmost_tnc_seq"]
    else:
        whichEnd = ["leftmost_tnc_seq"]

    BCDf = BaseCompAbCalc(SampleDf, group, whichEnd, lvl)

    # Save the results.
    BCDf.to_parquet("".join((args.output, "4_FrEIA/3_Abundances/",
                             group, "/", lvl, "/T__", filename, ".pq")),
                    compression="gzip",
                    index=False)


# Polynucleotide data is processed here.
def PolyNcAnalysis(SampleDf, args, filename, group, lvl):
    # Set ends needed to be analyzed.
    if args.doBoth:
        whichEnd = ["leftmost_read_seq", "rightmost_read_seq"]
    else:
        whichEnd = ["leftmost_read_seq"]

    OutDf = pd.DataFrame()
    h2PBCDf = pd. DataFrame()
    for e in whichEnd:
        # Convert string to characters in columns.
        hPBCDf = pd.DataFrame(SampleDf[e].apply(list).tolist())
        hPBCDf["Unit"] = SampleDf["chr"]
        hPBCDf = CastDataTypes(hPBCDf)
        #print(hPBCDf.info(memory_usage="deep"))
        h2PBCDf = pd.melt(hPBCDf,
                          id_vars=["Unit"],
                          value_vars=list(range(0, 10)),
                          var_name="Pos",
                          value_name="Base")
        h2PBCDf = CastDataTypes(h2PBCDf).reset_index(drop=True)
        h2PBCDf["Count"] = 1

        # Count bases in position.
        #print(h2PBCDf.info(memory_usage="deep"))
        OutDf = h2PBCDf.groupby(["Unit", "Pos", "Base"],
                                observed=True)["Count"].sum().reset_index()
        OutDf["CountSamp"] = (OutDf.groupby(["Pos", "Base"],
                                            observed=True)["Count"
                                                           ].transform(sum))
        OutDf["CountChr"] = (OutDf.groupby(["Unit", "Pos", "Base"],
                                           observed=True)["Count"
                                                          ].transform(sum))
        OutDf["RelAbSamp"] = (OutDf["CountSamp"]
                              / OutDf.groupby(["Pos"], observed=True
                                              )["Count"].transform(sum))
        OutDf["WhichEnd"] = e.split("_")[0]
        OutDf["WhichGroup"] = group
        OutDf.drop(["Count"], axis=1, inplace=True)

    # Save the results.
    OutDf.to_parquet("".join((args.output, "4_FrEIA/3_Abundances/",
                              group, "/", lvl, "/P__", filename, ".pq")),
                     compression="gzip", index=False)


# Joined ends data is processed here.
def JoinedEnds(SampleDf, args, filename, group, lvl):
    # Join mono- and trinucleotide ends.
    SampleDf["MJE"] = list(zip(SampleDf.leftmost_nc_seq,
                               SampleDf.rightmost_nc_seq))
    SampleDf["TJE"] = list(zip(SampleDf.leftmost_tnc_seq,
                               SampleDf.rightmost_tnc_seq))
    SampleDf["MJE"] = SampleDf["MJE"].astype("str")
    SampleDf["TJE"] = SampleDf["TJE"].astype("str")
    SampleDf["WhichEnd"] = "NA"

    # Cast string data into categorical for better memory usage
    # and faster calculations.
    for c in ["MJE", "TJE", "WhichEnd"]:
        SampleDf[c] = SampleDf[c].astype("category")
    JMEDf = BaseCompAbCalc(SampleDf, group, ["MJE"], lvl)
    JTEDf = BaseCompAbCalc(SampleDf, group, ["TJE"], lvl)

    # Save the results.
    JMEDf.to_parquet("".join((args.output, "4_FrEIA/3_Abundances/",
                              group, "/", lvl, "/JM__", filename, ".pq")),
                     compression="gzip",
                     index=False)
    JTEDf.to_parquet("".join((args.output, "4_FrEIA/3_Abundances/",
                              group, "/", lvl, "/JT__", filename, ".pq")),
                     compression="gzip",
                     index=False)


"""=== General ==="""


# Creating the output directories.
def CreateFolders(WhichGroup, whichLvl, args):
    for lvl in whichLvl:
        if not os.path.exists("".join((args.output, "4_FrEIA/3_Abundances/",
                                       WhichGroup, "/", lvl))):
            os.makedirs("".join((args.output, "4_FrEIA/3_Abundances/",
                                 WhichGroup, "/", lvl)))


# Loading the datasets from the extraction step


def ReadData(whichGroup, pathToFile, args):
    # print(pathToFile)
    rightEnds = ["rightmost_tnc_seq",
                 "rightmost_read_seq"]

    if pathToFile.split("/")[-1].endswith(".pq"):
        filename = pathToFile.split("/")[-1].strip(".pq")
        DataDf = pd.read_parquet(pathToFile)
    else:
        print("No datafiles found! Please supply the appropriate files!")

    # Cast string data into categorical for better memory usage
    # and faster calculations.
    # print(DataDf.info(memory_usage="deep"))
    DataDf = CastDataTypes(DataDf).reset_index(drop=True)
    # print(DataDf.info(memory_usage="deep"))

    DataDf.columns = Cols
    # Subsample dataset.
    if args.Subsample < 1:
        DataDf = DataDf.sample(frac=args.Subsample, replace=True)

    # Remove duplicated reads.
    DataDf.drop_duplicates(subset=["chr", "read_len",
                                   "leftmost_nc_pos",
                                   "rightmost_nc_pos"],
                           ignore_index=True,
                           inplace=True)

    # Remove rows containg N.
    DataDf = DataDf[~DataDf.leftmost_read_seq.str.contains("N")]
    DataDf = DataDf[~DataDf.rightmost_read_seq.str.contains("N")]

    # Keep reads with reasonable size.
    if args.fraLenMin:
        fraLenMin = args.fraLenMin
    else:
        fraLenMin = 1

    if args.fraLenMax:
        fraLenMax = args.fraLenMax
    else:
        fraLenMax = 500
    # print(DataDf[DataDf.read_len.isin(range(fraLenMin, fraLenMax))])
    DataDf = DataDf[DataDf.read_len.isin(range(fraLenMin, fraLenMax))]

    # Bin the datasets by genomic positions. This feature can be turned off
    # by 'subchroff' argument.
    if args.subchroff:
        DataDf = Binning(DataDf, args)
        whichLvl = ["Genome_Lvl", "Chromosome_Lvl", "SubChr_Lvl"]
    else:
        whichLvl = ["Genome_Lvl", "Chromosome_Lvl"]

    if args.SubSampRepeat and args.SubSampRepeat > 0:
        # Bootstrap sample the reads, and create N random tables.
        for n in range(0, args.SubSampRepeat):
            BSSampDf = DataDf.sample(frac=1, replace=True
                                     ).reset_index(drop=True)
            filename_bootstrap = "".join((filename, "_", str(n+1)))
            for lvl in whichLvl:
                MonoNcAnalysis(BSSampDf, args, filename_bootstrap,
                               whichGroup, lvl)
                TriNcAnalysis(BSSampDf, args, filename_bootstrap,
                              whichGroup, lvl)
                PolyNcAnalysis(BSSampDf, args, filename_bootstrap,
                               whichGroup, lvl)
                JoinedEnds(BSSampDf, args, filename_bootstrap,
                           whichGroup, lvl)
    else:
        for lvl in whichLvl:
            MonoNcAnalysis(DataDf, args, filename, whichGroup, lvl)
            TriNcAnalysis(DataDf, args, filename, whichGroup, lvl)
            # PolyNcAnalysis(DataDf, args, filename, whichGroup, lvl)
            JoinedEnds(DataDf, args, filename, whichGroup, lvl)


def Main():
    args = ParsingArguments()
    sampTDf = pd.read_csv(args.SampleTable, delim_whitespace=True)
    # Extract sample name from path.
    sampleName = args.inPath.split("/")[-1].strip(".pq")
    # Extract the group of the sample from sample sheet.
    WhichGroup = sampTDf.loc[sampTDf["sample_name"] == sampleName, "group"
                             ].item().strip("*")

    if args.subchroff:  # Select analysis type.
        whichLvl = ["Genome_Lvl", "Chromosome_Lvl", "SubChr_Lvl"]
    else:
        whichLvl = ["Genome_Lvl", "Chromosome_Lvl"]

    CreateFolders(WhichGroup, whichLvl, args)
    ReadData(WhichGroup, args.inPath, args)


if __name__ == "__main__":
    Main()
