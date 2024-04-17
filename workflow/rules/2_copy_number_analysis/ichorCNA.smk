import pandas as pd
import os
import glob

configfile: "../../../config/config.yaml" #Set config file.

Samplesheet = pd.read_csv(config["Samplesheet"], sep=' ') #Read sample sheet in a dataframe.

if config["Trimmer"] in ["bbduk", "cutadapt"]:
    ProjDirName = config["ProjName"] + "/trimmed"
    tmp_dir = config["TmpDir"] + "/" + ProjDirName #Set TEMPDIR.

elif config["Trimmer"] == "none":
    ProjDirName = config["ProjName"] + "/untrimmed"
    tmp_dir = config["TmpDir"] + "/" + ProjDirName #Set TEMPDIR.

else:
    print("\nMissing trimming option in config file!\n")
    exit()

# Get the dir name for the last created conda environment.
if config["conda_dir"]:
    conda_dir_path = config["conda_dir"]
else:
    try:
        fileL = glob.glob('.snakemake/conda/*')
        conda_dir_path = max(fileL, key=os.path.getctime)
    except:
        print("Run snakemake with --conda-create-envs-only first!")
        conda_dir_path = ""

localrules: all_ichorCNA

rule all_ichorCNA:
    input:
        expand(config["OutPath"] + "/" + ProjDirName + "/2_ichorCNA/{tumor}/{tumor}.cna.seg",
            tumor = Samplesheet["sample_name"]),
        expand(tmp_dir + "/2_read_depth/{sample}.bin{binSize}.wig",
            sample = Samplesheet["sample_name"],
            binSize = config["binSize"])

rule filter_and_index:
    input:
        config["OutPath"] + "/" + ProjDirName + "/1_mapping/{sample}.bam"
    output:
        bam = tmp_dir + "/2_filter_and_index/{sample}.bam",
        bai = tmp_dir + "/2_filter_and_index/{sample}.bam.bai"
    threads: config["ThreadNr"]
    resources: time = 10
    conda: "../../envs/copy_number_analysis_env.yaml"
    benchmark:
        config["OutPath"] + "/benchmark/" + ProjDirName + "/2_filter_and_index/{sample}.tsv"
    shell:
        """
        samtools view {input} \
        -h `# include header` \
        -q 5 `#> only reads with high quality` \
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

rule read_counter:
    input:
        bam = tmp_dir + "/2_filter_and_index/{sample}.bam",
        bai = tmp_dir + "/2_filter_and_index/{sample}.bam.bai"
    output:
        tmp_dir + "/2_read_depth/{sample}.bin{binSize}.wig"
    params:
        readCounter = config["readCounterScript"],
        binSize = config["binSize"],
        qual = "20",
        chrs = config["readDepth_chrs"] # changing this to chr1...., chr22 will cause PoN script not to work
    resources: time = 2
    conda: "../../envs/copy_number_analysis_env.yaml"
    benchmark:
        config["OutPath"] + "/benchmark/" + ProjDirName + "/2_read_depth/{sample}.bin{binSize}.tsv"
    # resources: mem_mb=4000 #> 4GB is overestimate by author
    log:
        config["OutPath"] + "/logs/" + ProjDirName + "/2_read_depth/{sample}.bin{binSize}.log"
    shell:
        """
        {params.readCounter} {input.bam} \
        -c {params.chrs} \
        -w {params.binSize} \
        -q {params.qual} > {output} 2> {log}
        """

rule ichorCNA:
    input:
        tum = tmp_dir + "/2_read_depth/{tumor}.bin" + str(config["binSize"]) + ".wig"
    output:
        #corrDepth = config["OutPath"] + "/" + ProjDirName + "/2_ichorCNA/{tumor}/{tumor}.correctedDepth.txt",
        #param = "config["OutPath"] + "/" + ProjDirName + "/2_ichorCNA/{tumor}/{tumor}.params.txt",
        cna = config["OutPath"] + "/" + ProjDirName + "/2_ichorCNA/{tumor}/{tumor}.cna.seg",
        #segTxt = config["OutPath"] + "/" + ProjDirName + "/2_ichorCNA/{tumor}/{tumor}.seg.txt",
        #seg = config["OutPath"] + "/" + ProjDirName + "/2_ichorCNA/{tumor}/{tumor}.seg",
        #rdata = config["OutPath"] + "/" + ProjDirName + "/2_ichorCNA/{tumor}/{tumor}.RData"
    params:
        outDir=config["OutPath"] + "/" + ProjDirName + "/2_ichorCNA/{tumor}/",
        rscript=conda_dir_path + config["ichorCNA_rscript"],
        id="{tumor}",
        ploidy=config["ichorCNA_ploidy"],
        normal=config["ichorCNA_normal"],
        gcwig=conda_dir_path + config["ichorCNA_gcWig"],
        mapwig=conda_dir_path + config["ichorCNA_mapWig"],
        normalpanel=config["ichorCNA_normalPanel"],
        estimateNormal=config["ichorCNA_estimateNormal"],
        estimatePloidy=config["ichorCNA_estimatePloidy"],
        scStates=config["ichorCNA_scStates"],
        maxCN=config["ichorCNA_maxCN"],
        includeHOMD=config["ichorCNA_includeHOMD"],
        chrs=config["ichorCNA_chrs"],
        chrTrain=config["ichorCNA_chrTrain"],
        genomeBuild=config["ichorCNA_genomeBuild"],
        genomeStyle=config["ichorCNA_genomeStyle"],
        centromere=conda_dir_path + config["ichorCNA_centromere"],
        fracReadsChrYMale=config["ichorCNA_fracReadsInChrYForMale"],
        minMapScore=config["ichorCNA_minMapScore"],
        exons=config["ichorCNA_exons"],
        txnE=config["ichorCNA_txnE"],
        txnStrength=config["ichorCNA_txnStrength"],
        plotFileType=config["ichorCNA_plotFileType"],
        plotYlim=config["ichorCNA_plotYlim"],
        libdir=config["ichorCNA_libdir"]
    resources:
        # mem_mb=10000 #> 10GB req is large overestimate by author!
        time = 5
    conda: "../../envs/copy_number_analysis_env.yaml"
    benchmark:
        config["OutPath"] + "/benchmark/" + ProjDirName + "/2_ichorCNA/{tumor}.tsv"
    log:
        config["OutPath"] + "/logs/" + ProjDirName + "/2_ichorCNA/{tumor}.log"
    shell:
        """
        Rscript {params.rscript} \
        --id {params.id} \
        --libdir {params.libdir} \
        --WIG {input.tum} \
        --gcWig {params.gcwig} \
        --mapWig {params.mapwig} \
        --normalPanel {params.normalpanel} \
        --ploidy \"{params.ploidy}\" \
        --normal \"{params.normal}\" \
        --maxCN {params.maxCN} \
        --includeHOMD {params.includeHOMD} \
        --chrs \"{params.chrs}\" \
        --chrTrain \"{params.chrTrain}\" \
        --genomeStyle {params.genomeStyle} \
        --genomeBuild {params.genomeBuild} \
        --estimateNormal {params.estimateNormal} \
        --estimatePloidy {params.estimatePloidy} \
        --estimateScPrevalence FALSE \
        --scStates \"{params.scStates}\" \
        --centromere {params.centromere} \
        --exons.bed {params.exons} \
        --txnE {params.txnE} \
        --txnStrength {params.txnStrength} \
        --minMapScore {params.minMapScore} \
        --fracReadsInChrYForMale {params.fracReadsChrYMale} \
        --plotFileType {params.plotFileType} \
        --plotYLim \"{params.plotYlim}\" \
        --outDir {params.outDir} > {log} 2> {log}
        """
