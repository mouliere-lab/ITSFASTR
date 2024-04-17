import pandas as pd

configfile: "../../../config/config.yaml"  # Set config file.

Samplesheet = pd.read_csv(config["Samplesheet"], sep='\s+') # Read sample sheet in a dataframe.

# Select output path depending on run mode.
metaPath = ""
print("Running in normal mode!")

if config["Trimmer"] in ["bbduk", "cutadapt"]:
    ProjDirName = config["ProjName"] + "/trimmed" + metaPath
    tmp_dir = config["TmpDir"] + "/" + ProjDirName + metaPath  # Set TEMPDIR.

elif config["Trimmer"] == "none":
    ProjDirName = config["ProjName"] + "/untrimmed" + metaPath
    tmp_dir = config["TmpDir"] + "/" + ProjDirName + metaPath  # Set TEMPDIR.

rule all_fragmentation_global:
        input:
            expand(config["OutPath"] + "/" + ProjDirName + "/3_fragmentation/{sample}_insert_size_metrics.txt",
                sample=Samplesheet["sample_name"])

if config["Mode"] == "ONT":  # Choose the sequencer used for data.
    rule ont_insert_size:
        input:
            config["OutPath"] + "/" + ProjDirName + "/1_mapping_quality/{sample}/NanoPlot-data.tsv.gz"
        output:
            tmp_dir + "/3_fragmentation/{sample}_insert_size_metrics.txt"
        shell:
           """
           gzip -dc {input} > {output}
           """

    rule ont_insert_size_simplify:
        input:
            tmp_dir + "/3_fragmentation/{sample}_insert_size_metrics.txt"
        output:
            config["OutPath"] + "/" + ProjDirName + "/3_fragmentation/{sample}_insert_size_metrics.txt"
        params:
            mapQ = config["FiltQual"]
        conda: "../../envs/fragmentation_env.yaml"
        shell:
            """
            python3 ../../scripts/3_fragmentation/insert_size_LR.py \
                    -i {input} \
                    -m {params.mapQ} \
                    -o {output}
            """

else:
    rule filter_and_index_fragmentation_global:
        input:
            config["OutPath"] + "/" + ProjDirName + "/1_mapping/{sample}.bam"
        output:
            bam = tmp_dir + "/2_filter_and_index/{sample}.bam",
            bai = tmp_dir + "/2_filter_and_index/{sample}.bam.bai"
        threads: config["ThreadNr"]
        params:
            qual = config["FiltQual"]
        conda: "../../envs/fragmentation_env.yaml"
        benchmark:
            config["OutPath"] + "/benchmark/" + ProjDirName + "/3_filter_and_index/{sample}.txt"
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

    rule insert_size:
        input:
            bam = tmp_dir + "/2_filter_and_index/{sample}.bam",
            bai = tmp_dir + "/2_filter_and_index/{sample}.bam.bai"  # Use the same bam and bai files that were previously created.
        output:
            txt = config["OutPath"] + "/" + ProjDirName + "/3_fragmentation/{sample}_insert_size_metrics.txt",
            pdf = tmp_dir+"/3_insert_size/{sample}_insert_size_histogram.pdf"
        conda: "../../envs/fragmentation_env.yaml"
        benchmark:
            config["OutPath"] + "/benchmark/" + ProjDirName + "/3_insert_size/{sample}.tsv"
        log:
            config["OutPath"] + "/logs/" + ProjDirName + "/3_insert_size/{sample}.log"
        shell:
            """
            picard CollectInsertSizeMetrics \
            INPUT={input.bam} \
            HISTOGRAM_WIDTH=1000 \
            MINIMUM_PCT=0.05 `# Default value` \
            OUTPUT={output.txt} \
            HISTOGRAM_FILE={output.pdf} `# Required parameter` \
            &> {log}
            """
