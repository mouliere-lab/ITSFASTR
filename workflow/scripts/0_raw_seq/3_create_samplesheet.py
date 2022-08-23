#!/usr/bin/env python3
# Norbert Moldován
import os
import argparse

#Parsing the arguments + some nice help text.
def ParsingArguments():
    parser=argparse.ArgumentParser()
    parser.add_argument("-d","--directory", dest = "Directory",
                        type=str, required=True, help="The path "
                        "to the folder to list the files from.")
    parser.add_argument("-o","--output", dest = "Output",
                        type=str, required=True, help="The path "
                        "to output samples sheet.")
    return parser.parse_args()

def Main():
    args=ParsingArguments()

    samplenameL = list()

    for root, dirs, files in os.walk(args.Directory):
        for filename in files:
            samplenameL.append(filename.split("_R")[0]) #Create a list of filenames.

    with open(args.Output, 'w') as filehandle:
        filehandle.write('%s\n' % "sample_name")
        for name in set(samplenameL):
            filehandle.write('%s\n' % name) #Print sample names to file.

    print("\n Samplesheet can be found at", args.Output + ".\n")

if __name__ == "__main__":
    Main()
