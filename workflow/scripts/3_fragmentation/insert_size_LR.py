#!/usr/bin/env python3
# Author: Norbert Moldován
# Mouliere Lab
# Amsterdam UMC
# Vrije Universiteit Amsterdam

import argparse
import pandas as pd

## This script transforms the output of rule ont_insert_size to a more managable smaller format
## similar to the output from PE Illumina data.

# Parsing the arguments + some nice help text.
def ParsingArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input",
                        dest="inPath",
                        type=str,
                        required=True,
                        help="Path to the incert size files.")
    parser.add_argument("-m", "--mapQ",
                        dest="mapQ",
                        type=int,
                        required=True,
                        help="The minimum mapping quality.")
    parser.add_argument("-o", "--output",
                        dest="output",
                        type=str,
                        required=True,
                        help="Output path.")
    return parser.parse_args()

def Main():
    args = ParsingArguments()
    insert_size_Df = pd.read_csv(args.inPath, delim_whitespace=True)
    insert_size_Df = insert_size_Df.loc[insert_size_Df["mapQ"] >= args.mapQ]
    simple_insert_size_Df = pd.DataFrame(insert_size_Df["aligned_lengths"].value_counts()).sort_index().reset_index()
    simple_insert_size_Df.rename(columns={"index": "insert_size",
                                          "aligned_lengths": "All_Reads.fr_count"},
                                  inplace=True)
    placeholders_Df = pd.DataFrame({"insert_size": ["### These lines are placeholders for cross-compatibility.",
                                                                                                  "###",
                                                                                                  "###",
                                                                                                  "###",
                                                                                                  "###",
                                                                                                  "###",
                                                                                                  "###",
                                                                                                  "###",
                                                                                                  "###",
                                                                                                  "insert_size"],
                                    "All_Reads.fr_count": ["###",
                                                           "###",
                                                           "###",
                                                           "###",
                                                           "###",
                                                           "###",
                                                           "###",
                                                           "###",
                                                           "###",
                                                           "All_Reads.fr_count"
                                        ]})
    simple_insert_size_Df = placeholders_Df.append(simple_insert_size_Df)
    simple_insert_size_Df.to_csv(args.output,
                                 index=False,
                                 sep="\t")

if __name__ == "__main__":
    Main()