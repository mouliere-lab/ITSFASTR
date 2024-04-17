import pandas as pd

configfile: "../../../config/config.yaml" #Set config file.

Samplesheet = pd.read_csv(config["Samplesheet"], sep=' ') #Read sample sheet in a dataframe.

localrules: all_mapping

# Select output path depending on run mode.
RefPath = config["RefPath"]
metaPath = ""
print("Running in normal mode!")

if config["Trimmer"] in ["bbduk", "cutadapt"]:
    ProjDirName = config["ProjName"] + "/trimmed" + metaPath
    tmp_dir = config["TmpDir"] + "/" + ProjDirName + metaPath  # Set TEMPDIR.
    MapIn = config["TmpDir"] + "/" + config["ProjName"] + "/trimmed" + "/1_trimming/"

elif config["Trimmer"] == "none":
    ProjDirName = config["ProjName"] + "/untrimmed" + metaPath
    tmp_dir = config["TmpDir"] + "/" + ProjDirName + metaPath  # Set TEMPDIR.
    MapIn = config["SamplePath"] + "/"

else:
    print("\nMissing trimming option in config file!\n")
    exit()

rule all_mapping:
    input:
        config["OutPath"] + "/" + ProjDirName + "/1_mapping_quality/multiqc/multiqc_report_interactive.html",
        expand(config["OutPath"] + "/" + ProjDirName + "/1_mapping/{sample}.bam",
                sample = Samplesheet["sample_name"])

if config["Mode"] == "ONT":  # Choose the sequencer used for data.
    rule map:
        input:
            MapIn + "{sample}.fastq.gz"
        output:
            tmp_dir + metaPath + "/1_mapping/{sample}.bam"
        params:
            ref = RefPath
        threads: config["ThreadNr"]
        conda: "../../envs/map_ONT_env.yaml"
        benchmark:
            config["OutPath"] + "/benchmark/" + ProjDirName + "/1_mapping/{sample}.tsv"
        log:
            config["OutPath"] + "/logs/" + ProjDirName + "/1_mapping/{sample}.log"
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
    rule map:
            input:
                MapIn + "{sample}_R1.fq.gz",
                MapIn + "{sample}_R2.fq.gz"
            output:
                tmp_dir + metaPath +"/1_mapping/{sample}.bam"
            params:
                ref = RefPath
            threads: config["ThreadNr"]
            resources: mem_mb=36000+8000
                       #> Prevent memory exhaustion by samtools sort
                       #>> 1500MB/thread for sorting, 8000 for BWA
            conda: "../../envs/preprocessing_env.yaml"
            benchmark:
                config["OutPath"] + "/benchmark/" + ProjDirName + "/1_mapping/{sample}.tsv"
            log:
                config["OutPath"] + "/logs/" + ProjDirName + "/1_mapping/{sample}.log"
            shell:
                """
                bwa mem {params.ref} {input} \
                -t {threads} \
                2> {log} | \
                samtools sort \
                -@ {threads} \
                -o {output} 2>> {log}
                """

if config["DupMarker"] == "picard_MD":
    rule mark_duplicates:
        input:
            tmp_dir + "/1_mapping/{sample}.bam"
        output:
            bam = config["OutPath"] + "/" + ProjDirName + "/1_mapping/{sample}.bam",
            mtrx = config["OutPath"] + "/logs/" + ProjDirName + "/2_mark_duplicates/{sample}_dup_metrics.txt"
        priority: 100
        threads: config["ThreadNr"]
        resources: max_file_handles=50000,
                   mem_mb=54000
        params: jar = "$CONDA_PREFIX/share/picard-2.22.2-0/picard.jar",
                TMP_DIR = tmp_dir + "/2_mark_duplicates/"
        conda: "../../envs/preprocessing_env.yaml"
        benchmark:
            config["OutPath"] + "/benchmark/" + ProjDirName + "/2_mark_duplicates/{sample}.tsv"
        log:
            config["OutPath"] + "/logs/" + ProjDirName + "/2_mark_duplicates/{sample}.log"
        shell:
            """
            java -Xmx{resources.mem_mb}M\
            -XX:ParallelGCThreads={threads}\
            -jar {params.jar} MarkDuplicates\
            INPUT={input}\
            REMOVE_DUPLICATES=false\
            MAX_FILE_HANDLES={resources.max_file_handles}\
            CREATE_INDEX=false\
            COMPRESSION_LEVEL=9\
            TMP_DIR={params.TMP_DIR}\
            METRICS_FILE={output.mtrx}\
            OUTPUT={output.bam}\
            &> {log}
            """
else:
    if config["DupMarker"] != "sambamba_MD":
        print("WARNING! The selected duplicate marker is not available! \n",
              "Running Sambamba markdup as default!")

    rule mark_duplicates:
        input:
            tmp_dir + metaPath + "/1_mapping/{sample}.bam"
        output:
            bam = config["OutPath"] + "/" + ProjDirName + "/1_mapping/{sample}.bam",
        threads: config["ThreadNr"]
        params:
            tmpDir = tmp_dir + "/2_mark_duplicates"
        conda: "../../envs/preprocessing_env.yaml"
        benchmark:
            config["OutPath"] + "/benchmark/" + ProjDirName + "/2_mark_duplicates/{sample}.tsv"
        log:
            config["OutPath"] + "/logs/" + ProjDirName + "/2_mark_duplicates/{sample}.log"
        shell:
            """
            sambamba markdup \
            -t {threads} \
            --tmpdir {params.tmpDir} \
            {input} \
            {output} \
            2>> {log}
            """

if config["MethLib"]:  # For converted libraries run biscuit QC.
    rule biscuit_QC:
        input:
            config["OutPath"] + "/" + ProjDirName + "/1_mapping/{sample}.bam"
        output:
            config["OutPath"] + "/" + ProjDirName + "/1_mapping_quality/biscuit/{sample}/biscuit_qc"
        params:
            ref = RefPath + ".fna"
        conda: "../../envs/biscuit_env.yaml"
        benchmark:
            config["OutPath"] + "/benchmark/" + ProjDirName + "/1_mapping_quality/biscuit/{sample}.tsv"
        shell:
            """
            biscuit qc {params.ref} {input} {output}
            """

rule flagstat:
    input:
        config["OutPath"] + "/" + ProjDirName + "/1_mapping/{sample}.bam"
    output:
        config["OutPath"] + "/" + ProjDirName + "/1_mapping_flagstat/{sample}_flagstats.txt"
    threads: config["ThreadNr"]
    resources:
        time = 10
    conda: "../../envs/preprocessing_env.yaml"
    benchmark:
        config["OutPath"] + "/benchmark/" + ProjDirName + "/1_mapping_flagstat/{sample}.tsv"
    shell:
        """
        samtools flagstat {input} --threads {threads} > {output}
        """

if config["Mode"] == "ONT":  # Choose the sequencer used for data.
    rule NanoPlot:
        input:
            config["OutPath"] + "/" + ProjDirName + "/1_mapping/{sample}.bam"
        output:
            directory(config["OutPath"] + "/" + ProjDirName + "/1_mapping_quality/{sample}")
        threads: config["ThreadNr"]
        conda: "../../envs/map_ONT_env.yaml"
        benchmark:
            config["OutPath"] + "/benchmark/" + ProjDirName + "/1_mapping_quality/{sample}.tsv"
        log:
            config["OutPath"] + "/logs/" + ProjDirName + "/1_mapping_quality/{sample}.log"
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
else:
    rule qualimap:
        input:
            config["OutPath"] + "/" + ProjDirName + "/1_mapping/{sample}.bam"
        output:
            dir=directory(config["OutPath"] + "/" + ProjDirName + "/1_mapping_quality/{sample}")
        threads: config["ThreadNr"]
        resources:
                mem_mb=60000, #round(60/(math.sqrt(12))),
                time=120
        conda: "../../envs/preprocessing_env.yaml"
        benchmark:
            config["OutPath"] + "/benchmark/" + ProjDirName + "/1_mapping_quality/{sample}.tsv"
        log:
            config["OutPath"] + "/logs/" + ProjDirName + "/1_mapping_quality/{sample}.log"
        shell:
            """
            unset DISPLAY; \
            qualimap bamqc -bam {input} \
            --java-mem-size={resources.mem_mb}M \
            -nt {threads} \
            -outdir {output.dir} \
            `# -outfile {output} #> gives errors` \
            &> {log}
            """

rule multiqc_mapping:
    input:
        expand(config["OutPath"] + "/" + ProjDirName + "/1_mapping_quality/{sample}",
        sample = Samplesheet["sample_name"]),
        expand(config["OutPath"] + "/" + ProjDirName + "/1_mapping_flagstat/{sample}_flagstats.txt",
        sample = Samplesheet["sample_name"]),
        expand(config["OutPath"] + "/logs/" + ProjDirName + "/2_mark_duplicates/{sample}.log",
        sample = Samplesheet["sample_name"])
    output:
        interactive = config["OutPath"] + "/" + ProjDirName + "/1_mapping_quality/multiqc/multiqc_report_interactive.html",
        flat = config["OutPath"] + "/" + ProjDirName + "/1_mapping_quality/multiqc/multiqc_report_flat.html"
    resources:
        time = 60
    conda: "../../envs/preprocessing_env.yaml"
    benchmark:
        config["OutPath"] + "/benchmark/" + ProjDirName + "/1_mapping_multiqc/multiqc.tsv"
    log:
        interactive = config["OutPath"] + "/logs/" + ProjDirName + "/1_mapping_multiqc/multiqc_interactive.log",
        flat = config["OutPath"] + "/logs/" + ProjDirName + "/1_mapping_multiqc/multiqc_flat.log"
    shell:
        """
        multiqc {input} \
        --interactive \
        --force \
        -o $(dirname {output.interactive}) \
        --filename $(basename {output.interactive}) `# Next-level MultiQC hack` \
        --verbose &> {log.interactive} &

        multiqc {input} \
        --flat \
        --force \
        -o $(dirname {output.flat}) \
        --filename $(basename {output.flat}) `# Next-level MultiQC hack` \
        --verbose &> {log.flat} &
        wait
        """
