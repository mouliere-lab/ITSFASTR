#!/usr/bin/env python3
# Norbert Moldován
import os
import argparse
import subprocess

#Parsing the arguments + some nice help text.
def ParsingArguments():
    parser=argparse.ArgumentParser()
    parser.add_argument("-d","--directory", dest = "Directory",
                        type=str, required=True, help="The path "
                        "to the folder with the files.")
    parser.add_argument("-l","--lane_indicator", dest = "LaneIndicator",
                        type=str, required=True, help="The pattern "
                        "indicating the lanes, separated by a coma. Ex: _L001,_L002")
    parser.add_argument("-o","--output", dest = "Output",
                        type=str, required=True, help="The path "
                        "to output the merged files to.")
    return parser.parse_args()

def Main():
    args=ParsingArguments()

    samplenameS = set()

    for root, dirs, files in os.walk(args.Directory):
        for filename in files:
            for lane in args.LaneIndicator.split(","):
                if lane in filename:
                    samplenameS.add(filename.replace(lane,"*")) # Create a list of filenames without the lane indicator string.

    for samplename in samplenameS:
        cat_fq = 'cat {} > {}'.format("".join((args.Directory, "/", samplename)), "".join((args.Output, "/", samplename.replace("*",""))))
        subprocess.run(cat_fq, shell=True)

if __name__ == "__main__":
    Main()
