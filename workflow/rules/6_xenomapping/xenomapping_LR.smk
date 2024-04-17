import pandas as pd

configfile: "../../../config/config.yaml"  # Set config file.

Samplesheet = pd.read_csv(config["Samplesheet"], sep=' ') # Read sample sheet in a dataframe.

# Select output path depending on run mode.
metaPath = ""
print("Running in normal mode!")

if config["Trimmer"] in ["bbduk", "cutadapt"]:
    ProjDirName = config["ProjName"] + "/trimmed" + metaPath
    tmp_dir = config["TmpDir"] + "/" + ProjDirName + metaPath  # Set TEMPDIR.

elif config["Trimmer"] == "none":
    ProjDirName = config["ProjName"] + "/untrimmed" + metaPath
    tmp_dir = config["TmpDir"] + "/" + ProjDirName + metaPath  # Set TEMPDIR.


rule xenomapping_all:
    input:
        expand(config["OutPath"] + "/" + ProjDirName + "/1_mapping/{sample}_{ps}.bam",
               sample=Samplesheet["sample_name"],
               ps=["primary",
                   "secondary"]),
        expand(config["OutPath"] + "/" + ProjDirName + "/1_mapping_quality/{sample}_{ps}",
                         sample=Samplesheet["sample_name"],
                         ps=["primary",
                             "secondary"])

if config["Mode"] == "ONT":  # Choose the sequencer used for data.
    rule map_primary:
        input:
            tmp_dir + metaPath + "/1_trimming/{sample}.fastq.gz"
        output:
            tmp_dir + metaPath + "/7_xenomapping/1_mapping/{sample}_primary.bam"
        params:
            ref = config["RefPath_Primary"]
        threads: config["ThreadNr"]
        conda: "../../envs/map_ONT_env.yaml"
        benchmark:
            config["OutPath"] + "/benchmark/" + ProjDirName + "/7_xenomapping/1_mapping/{sample}_primary.tsv"
        log:
            config["OutPath"] + "/logs/" + ProjDirName + "/7_xenomapping/1_mapping/{sample}_primary.log"
        shell:
            """
            minimap2 -ax map-ont {params.ref} \
                     {input} \
                     -t {threads} \
                     2> {log} |
                     samtools sort \
                     -@ {threads} \
                     -o {output} 2>> {log}
            """
            
    
    rule map_secondary:
        input:
            tmp_dir + metaPath + "/1_trimming/{sample}.fastq.gz"
        output:
            tmp_dir + metaPath + "/7_xenomapping/1_mapping/{sample}_secondary.bam"
        params:
            ref = config["RefPath_Secondary"]
        threads: config["ThreadNr"]
        conda: "../../envs/map_ONT_env.yaml"
        benchmark:
            config["OutPath"] + "/benchmark/" + ProjDirName + "/7_xenomapping/1_mapping/{sample}_secondary.tsv"
        log:
            config["OutPath"] + "/logs/" + ProjDirName + "/7_xenomapping/1_mapping/{sample}_secondary.log"
        shell:
            """
            minimap2 -ax map-ont {params.ref} \
                     {input} \
                     -t {threads} \
                     2> {log} |
                     samtools sort \
                     -@ {threads} \
                     -o {output} 2>> {log}
            """

else:
    print("This script is used to separate long ONT or PacBio reads. To separate Illumina reads please use another tool!")
    exit()


rule filter_read_names:
    input:
        tmp_dir + metaPath + "/7_xenomapping/1_mapping/{sample}_{ps}.bam"
    output:
        tmp_dir + metaPath + "/7_xenomapping/filtered_names/{sample}_{ps}.txt"
    threads: config["ThreadNr"]
    params:
        map_score = config["min_score"]
    conda: "../../envs/map_ONT_env.yaml"
    benchmark:
        config["OutPath"] + "/benchmark/" + ProjDirName + "/7_xenomapping/filtered_names/{sample}_{ps}.tsv"
    log:
        config["OutPath"] + "/logs/" + ProjDirName + "/7_xenomapping/filtered_names/{sample}_{ps}.log"
    shell:
        """
        samtools view {input} \
        -q {params.map_score} `#> only reads with high quality` \
        -F 4 `# not unmapped` \
        -@ {threads} |\
        cut -f 1,5 `# select qnames and mapq`|\
        sort -k1,1 | \
        uniq -f 1 \
        > {output} 2>> {log}
        """


# Thise rules select alignments unique to the primary and secondary files.
rule unique_names:
    input:
        namesP = tmp_dir + metaPath + "/7_xenomapping/filtered_names/{sample}_primary.txt",
        namesS = tmp_dir + metaPath + "/7_xenomapping/filtered_names/{sample}_secondary.txt"
    output:
        namesP_U = tmp_dir + metaPath + "/7_xenomapping/unique_names/{sample}_primary.txt",
        namesS_U = tmp_dir + metaPath + "/7_xenomapping/unique_names/{sample}_secondary.txt"
    benchmark:
        config["OutPath"] + "/benchmark/" + ProjDirName + "/7_xenomapping/unique_names/{sample}_primary.tsv"
    log:
        config["OutPath"] + "/logs/" + ProjDirName + "/7_xenomapping/unique_names/{sample}_primary.log"
    shell:
        """
        awk 'BEGIN {{OFS="\t"}} NR==FNR {{a[$1]=$2; next}} {{if ($1 in a) {{if ($2>a[$1]) {{delete a[$1]}}}} }} END {{for (i in a) print i}}' {input.namesP} {input.namesS} \
        > {output.namesP_U} \
        && \
        awk 'BEGIN {{OFS="\t"}} NR==FNR {{a[$1]=$2; next}} {{if ($1 in a) {{if ($2>a[$1]) {{delete a[$1]}}}} }} END {{for (i in a) print i}}' {input.namesS} {input.namesP} \
        > {output.namesS_U} \
        2>> {log}
        """


rule select_alignments_primary:
    input:
        bamP = tmp_dir + metaPath + "/7_xenomapping/1_mapping/{sample}_primary.bam",
        namesP_U = tmp_dir + metaPath + "/7_xenomapping/unique_names/{sample}_primary.txt",
        namesS_U = tmp_dir + metaPath + "/7_xenomapping/unique_names/{sample}_secondary.txt"
    output:
        bamP_U = tmp_dir + metaPath + "/7_xenomapping/unique_mapping/{sample}_primary.bam"
    threads: config["ThreadNr"]
    conda: "../../envs/map_ONT_env.yaml"
    benchmark:
        config["OutPath"] + "/benchmark/" + ProjDirName + "/7_xenomapping/unique_names/{sample}_primary.tsv"
    log:
        logP = config["OutPath"] + "/logs/" + ProjDirName + "/7_xenomapping/unique_names/{sample}_primary.log"
    shell:
        """
        samtools view {input.bamP} \
        -h `# include header` \
        -N {input.namesP_U} \
        -@ {threads} \
        -b |\
        samtools sort -\
        -@ {threads} \
        -o {output.bamP_U} 2>> {log.logP} \
        """


rule select_alignments_secondary:
    input:
        bamS = tmp_dir + metaPath + "/7_xenomapping/1_mapping/{sample}_secondary.bam",
        namesP_U = tmp_dir + metaPath + "/7_xenomapping/unique_names/{sample}_primary.txt",
        namesS_U = tmp_dir + metaPath + "/7_xenomapping/unique_names/{sample}_secondary.txt"
    output:
        bamS_U = tmp_dir + metaPath + "/7_xenomapping/unique_mapping/{sample}_secondary.bam"
    threads: config["ThreadNr"]
    conda: "../../envs/map_ONT_env.yaml"
    benchmark:
        config["OutPath"] + "/benchmark/" + ProjDirName + "/7_xenomapping/unique_names/{sample}_primary.tsv"
    log:
        logS = config["OutPath"] + "/logs/" + ProjDirName + "/7_xenomapping/unique_names/{sample}_secondary.log"
    shell:
        """
        samtools view {input.bamS} \
        -h `# include header` \
        -N {input.namesS_U} \
        -@ {threads} \
        -b |\
        samtools sort -\
        -@ {threads} \
        -o {output.bamS_U} 2>> {log.logS}
        """


rule mark_duplicates:
    input:
        tmp_dir + metaPath + "/7_xenomapping/unique_mapping/{sample}_{ps}.bam"
    output:
        bam = config["OutPath"] + "/" + ProjDirName + "/1_mapping/{sample}_{ps}.bam",
    threads: config["ThreadNr"]
    params:
        tmpDir = tmp_dir + "/7_xenomapping/3_mark_duplicates/"
    conda: "../../envs/preprocessing_env.yaml"
    benchmark:
        config["OutPath"] + "/benchmark/" + ProjDirName + "/7_xenomapping/4_mark_duplicates/{sample}_{ps}.tsv"
    log:
        config["OutPath"] + "/logs/" + ProjDirName + "/7_xenomapping/4_mark_duplicates/{sample}_{ps}.log"
    shell:
        """
        sambamba markdup \
        -t {threads} \
        --tmpdir {params.tmpDir} \
        {input} \
        {output} \
        2>> {log}
        """


rule NanoPlot:
    input:
        config["OutPath"] + "/" + ProjDirName + "/1_mapping/{sample}_{ps}.bam"
    output:
        directory(config["OutPath"] + "/" + ProjDirName + "/1_mapping_quality/{sample}_{ps}")
    threads: config["ThreadNr"]
    conda: "../../envs/map_ONT_env.yaml"
    benchmark:
        config["OutPath"] + "/benchmark/" + ProjDirName + "/1_mapping_quality/{sample}_{ps}.tsv"
    log:
        config["OutPath"] + "/logs/" + ProjDirName + "/1_mapping_quality/{sample}_{ps}.log"
    shell:
        """
        NanoPlot --bam {input} \
                    -o {output} \
                    --raw \
                    --alength \
                    -t {threads} \
                    --huge \
                    2>> {log}
        """