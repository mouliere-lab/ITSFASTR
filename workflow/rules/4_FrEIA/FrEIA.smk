import pandas as pd
import numpy as np

configfile: "../../../config/config.yaml"  # Set config file.
# Read sample sheet in a dataframe.
Samplesheet = pd.read_csv(config["Samplesheet"], sep=' ')

ControlGroup = [c for c in Samplesheet["group"].unique() if "*" in c]
AffectedGroups = ",".join(np.setdiff1d(Samplesheet["group"].unique(),
                                       ControlGroup))
ControlGroup = ControlGroup[0][:-1]
SampleGroups = Samplesheet["group"].str.replace(r'\W', '')

Levels = ["Genome_Lvl"]  # Only genome and chromosome lvl analysis is done.
Prefixes = ["M", "T"]

metaPath = ""

# Collecting the outut paths for the multioutput rule FrEIA_CalcAb.
def CalcAbOut(group, sample):
    outPathL = list()
    for lvl in Levels:
        for pref in Prefixes:
            outPath = "".join((config["OutPath"], "/", ProjDirName, metaPath,
                               "/4_FrEIA/3_Abundances/", group, "/",
                               lvl, "/", pref, "__", sample, ".pq"))
            outPathL.append(outPath)
    return outPathL

if config["Trimmer"] in ["bbduk", "cutadapt"]:
    ProjDirName = config["ProjName"] + "/trimmed" + metaPath
    tmp_dir = config["TmpDir"] + "/" + ProjDirName + metaPath  # Set TEMPDIR.

elif config["Trimmer"] == "none":
    ProjDirName = config["ProjName"] + "/untrimmed" + metaPath
    tmp_dir = config["TmpDir"] + "/" + ProjDirName + metaPath  # Set TEMPDIR.

localrules: compare_groups

rule compare_groups:
    input:
        expand(config["OutPath"] + "/" + ProjDirName +
               "/4_FrEIA/3_Abundances/{group}/{lvl}/{prefix}__{sample}.pq",
               zip,
               group=np.repeat(SampleGroups, (len(Prefixes) * len(Levels))),
               lvl=np.repeat(Levels, len(Prefixes)).tolist() * len(SampleGroups),
               prefix=Prefixes * len(Levels) * len(SampleGroups),
               sample=np.repeat(Samplesheet["sample_name"],
                                (len(Prefixes) * len(Levels))).tolist() * len(SampleGroups))
    threads: config["ThreadNr"]
    params:
        outPath = config["OutPath"] + "/" + ProjDirName + "/",
        sampTable = config["Samplesheet"],
        subsResults = config["SubsetResults"],
        regroup = config["Regroup"]
    conda: "../../envs/FrEIA_env.yaml"
    shell:
        """
        python3 ../../scripts/4_FrEIA/4_compare_groups.py \
        -i {params.outPath} \
        -st {params.sampTable} \
        -t {threads} \
        -sN {params.subsResults} \
        -rgr {params.regroup} \
        #--plot
        """

# Calculating base and motif fractions per sample.
rule data_transformation:
    input:
        (config["OutPath"] + "/" + ProjDirName +
         "/4_FrEIA/1_extract_fragment_ends/{sample}.pq")
    output:
        CalcAbOut("{group}", "{sample}")
    params:
        outPath = config["OutPath"] + "/" + ProjDirName + "/",
        sampTable = config["Samplesheet"],
        fraSizeMin = config["FragmSizeMin"],
        fraSizeMax = config["FragmSizeMax"],
        subSamp = config["SubSampleRate"],
        bsSampNr = config["BsSampNr"]
    conda: "../../envs/FrEIA_env.yaml"
    log:
        (config["OutPath"] + "/logs/" + ProjDirName +
         "/4_FrEIA/3_data_transformation/{group}/{sample}.log")
    shell:
        """
        python3 ../../scripts/4_FrEIA/3_data_transformation.py \
        -i {input} \
        -o {params.outPath} \
        -st {params.sampTable} \
        -fsmin {params.fraSizeMin} \
        -fsmax {params.fraSizeMax} \
        -subs {params.subSamp} \
        --bootstrap_sample {params.bsSampNr} \
        --subchroff \
        2> {log}
        """
