#!/usr/bin/env python3
# Author: Norbert Moldován

import os, sys, re
import argparse
import time
from multiprocessing import Pool
from functools import partial
import pandas as pd
import pysam
from Bio.SeqUtils import GC
from FrEIA_tools import CastDataTypes


def ParsingArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--bam",
                        dest="bamIn",
                        type=str,
                        required=True,
                        help="Path to the input bam file.")
    parser.add_argument("-n", "--nrBase",
                        dest="nrBase",
                        type=int,
                        required=False,
                        help="Number of bases to fetch. Default is 10.")
    parser.add_argument("-t", "--threads",
                        dest="threads",
                        type=int,
                        required=True,
                        help="Number of threads.")
    parser.add_argument("-o", "--out",
                        dest="output",
                        type=str,
                        required=False,
                        help="Path to the output file.")
    parser.add_argument("-c", "--contigs",
                        dest="contigs",
                        type=str,
                        required=True,
                        help="List of contig names separated by a coma.")
    parser.add_argument("-m", "--mode",
                        dest="mode",
                        type=str,
                        required=True,
                        help="Chose the sequencer [Illumina, ONT]."
                             "Default: Illumina.")
    return parser.parse_args()


def getReadSeq(read, mate, nrBase):

    if not read.is_reverse:
        P1_seq = read.query_alignment_sequence[:nrBase]
        P2_seq = mate.query_alignment_sequence[-nrBase:][::-1]  # Rev sequence.
    else:
        P2_seq = read.query_alignment_sequence[-nrBase:][::-1]  # Rev sequence.
        P1_seq = mate.query_alignment_sequence[:nrBase]

    GC_cont = GC(read.query_alignment_sequence + mate.query_alignment_sequence)

    return {"P1_seq": P1_seq,
            "P2_seq": P2_seq,
            "GC_cont": GC_cont}


def seLength(cigarStr):
    cig_dict = {"M" : 0, "I" : 0, "D" : 0, "N" : 0}
    cigar_list = []
    for number, event, in re.findall("(\d+)([ISMDHN])", cigarStr):
        cigar_list.append((event, number),)
        if event == "M":
            cig_dict["M"] += int(number)  # Mapped count.
        elif event == "I":
            cig_dict["I"] += int(number)  # Insert count.
        elif event == "D":
            cig_dict["D"] += int(number)  # Del count.
        elif event == "N":
            cig_dict["N"] += int(number)

    q_length = int(cig_dict["M"])+int(cig_dict["D"])-1
    return q_length


def fetchData(args, chromosome):

    if args.nrBase:
        nrBase = args.nrBase
    else:
        nrBase = 10  # Defoult value for nrBase.

    outDict = {"qname": [],
               "chr": [],
               "read_len": [],
               "strand": [],
               "leftmost_nc_pos": [],
               "rightmost_nc_pos": [],
               "leftmost_nc_seq": [],
               "leftmost_tnc_seq": [],
               "rightmost_tnc_seq": [],
               "rightmost_nc_seq": [],
               "leftmost_read_seq": [],
               "rightmost_read_seq": [],
               "gccont": []}

    bamfile = pysam.AlignmentFile(args.bamIn, "rb")
    readCache = {}

    # readCount = bamfile.count(chromosome)
    # i = 0
    for read in bamfile.fetch(chromosome):
        # print(chromosome, ":", ((i)/readCount)*100,
        #        "Cache size:", len(readCache), end='\r')
        # i += 1
        if args.mode == "ONT":  # Delete cache every round for ONT seq.
            readCache = {}
            readCache[read.query_name] = read
        elif len(readCache) > 10000:  # Limit cache size if not ONT seq.
            readCache = {}
        # Cache reads as Pysam's mate method is too slow.
        if read.query_name not in readCache:
            readCache[read.query_name] = read
            continue
        else:
            mate = readCache[read.query_name]
            SeqDict = getReadSeq(read, mate, nrBase)  # Get the sequences.

            if ((SeqDict.get("P1_seq") != "") &
               (SeqDict.get("P2_seq") != "")):
                outDict["qname"].append(read.query_name)
                outDict["chr"].append(bamfile.get_reference_name(read.reference_id))
                if args.mode == "ONT":
                    outDict["read_len"].append(seLength(read.cigarstring))
                else:
                    outDict["read_len"].append(abs(read.template_length))
                outDict["strand"].append(0)
                outDict["leftmost_nc_pos"].append(read.reference_start)
                outDict["rightmost_nc_pos"].append(read.reference_start
                                                   + read.template_length)

                outDict["leftmost_nc_seq"].append(SeqDict.get("P1_seq")[0])
                outDict["rightmost_nc_seq"].append(SeqDict.get("P2_seq")[0])
                outDict["leftmost_tnc_seq"].append(SeqDict.get("P1_seq")[:3])
                outDict["rightmost_tnc_seq"].append(SeqDict.get("P2_seq")[:3])
                outDict["leftmost_read_seq"].append(SeqDict.get("P1_seq"))
                outDict["rightmost_read_seq"].append(SeqDict.get("P2_seq"))
                outDict["gccont"].append(SeqDict.get("GC_cont"))

        del readCache[read.query_name]

    bamfile.close()

    if len(pd.DataFrame.from_dict(outDict)) > 0:
        return CastDataTypes(pd.DataFrame.from_dict(outDict))
    else:
        return pd.DataFrame.from_dict(outDict)


def Main():
    # Tt = time.time()
    args = ParsingArguments()

    ChrL = args.contigs.split(",")

    OutDf = pd.DataFrame()

    # print("File size:", os.stat(args.bamIn).st_size / (1024*1024), "MB")
    # For larg inputs threading is turned off to conserv memory.
    if (os.stat(args.bamIn).st_size / (1024*1024)) > 8000:
        args.threads = 1

    pool = Pool(processes=args.threads)

    OutDf = OutDf.append(pool.map(partial(fetchData, args), ChrL),
                         ignore_index=True)

    pool.close()
    pool.join()

    #print(OutDf.info(memory_usage="deep"))

    #print(OutDf["leftmost_nc_seq"].value_counts())
    OutDf.to_parquet(args.output,  # Save the results.
                     compression="gzip",
                     index=False)
    # print(args.bamIn,": ",time.time()-Tt, "\n")


if __name__ == "__main__":
    Main()
