import pandas as pd

configfile: "../../../config/config.yaml"  # Set config file.
# Read sample sheet in a dataframe.
Samplesheet = pd.read_csv(config["Samplesheet"], sep=' ')

# Select output path depending on run mode.
metaPath = ""
if config["Trimmer"] in ["bbduk", "cutadapt"]:
    ProjDirName = config["ProjName"] + "/trimmed" + metaPath
    tmp_dir = config["TmpDir"] + "/" + ProjDirName + metaPath  # Set TEMPDIR.

elif config["Trimmer"] == "none":
    ProjDirName = config["ProjName"] + "/untrimmed" + metaPath
    tmp_dir = config["TmpDir"] + "/" + ProjDirName + metaPath  # Set TEMPDIR.

localrules: all_preprocessing

report: (config["OutPath"] + "/" + ProjDirName + "reports/4_FrEIA_preprocessing.rst")

rule all_preprocessing:
    input:
        expand(config["OutPath"] + "/" + ProjDirName +
               "/4_FrEIA/1_extract_fragment_ends/{sample}.pq",
               sample=Samplesheet["sample_name"])

rule filter_and_index_FrEIA:
    input:
        config["OutPath"] + "/" + ProjDirName + "/1_mapping/{sample}.bam"
    output:
        bam = tmp_dir + "/2_filter_and_index/{sample}.bam",
        bai = tmp_dir + "/2_filter_and_index/{sample}.bam.bai"
    threads: config["ThreadNr"]
    params:
        qual = config["FiltQual"]
    conda: "../../envs/FrEIA_env.yaml"
    benchmark:
        (config["OutPath"] + "/benchmark/" +
         ProjDirName + "/2_filter_and_index/{sample}.tsv")
    shell:
        """
        samtools view {input} \
        -h `# include header` \
        -q {params.qual} `#> only reads with high quality` \
        -F 4 `# not unmapped` \
        -F 256 `# no secondary alignment; thus primary alignment` \
        -F 1024 `# no PCR duplicate` \
        -F 2048 `# no supplementary (=chimeric) alignment` \
        -@ {threads} \
        -b | tee {output.bam} | \
        \
        samtools index - `# stdin` \
        -@ {threads} \
        {output.bai}
        """

rule FragmEndSeq:
    input:
        tmp_dir + "/2_filter_and_index/{sample}.bam"
    output:
        report(config["OutPath"] + "/" + ProjDirName +
               "/4_FrEIA/1_extract_fragment_ends/{sample}.pq",
               category="Fragment end retrieval")
    threads: config["ThreadNr"]
    params:
        nBase = config["nBase"],
        contigs = config["Contigs"],
        mode = config["Mode"]
    conda: "../../envs/FrEIA_env.yaml"
    benchmark:
        (config["OutPath"] + "/benchmark/" + ProjDirName +
         "/4_FrEIA/1_extract_fragment_ends/{sample}.tsv")
    log:
        (config["OutPath"] + "/logs/" + ProjDirName +
         "/4_FrEIA/1_extract_fragment_ends/{sample}.log")
    shell:
        """
        python3 ../../scripts/4_FrEIA/1_extract_fragment_ends.py \
        -b {input} \
        -o {output} \
        -c {params.contigs} \
        -t {threads} \
        -n {params.nBase} \
        -m {params.mode} \
        2> {log}
        """
