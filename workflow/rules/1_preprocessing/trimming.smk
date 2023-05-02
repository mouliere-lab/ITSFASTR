import pandas as pd

configfile: "../../../config/config.yaml"  # Set config file.

ProjDirName = config["ProjName"] + "/trimmed"
tmp_dir = config["TmpDir"] + "/" + ProjDirName  # Set TEMPDIR.
# Read sample sheet in a dataframe.
Samplesheet = pd.read_csv(config["Samplesheet"], delim_whitespace=True)

localrules: all_trim_ONT, all_trim_Illumina

if config["Mode"] == "ONT":  # Choose the sequencer used for data.
    rule all_trim_ONT:
        input:
            expand(tmp_dir + "/1_trimming/{porechop_sample}.fastq.gz",
                   zip,
                   porechop_sample=Samplesheet["sample_name"])

    rule trim_porechop:
        input:
            config["SamplePath"] + "/{porechop_sample}"
        output:
            tmp_dir + "/1_trimming/{porechop_sample}.fastq.gz"
        threads: config["ThreadNr"]
        conda: "../../envs/porechop_env.yaml"
        benchmark:
            (config["OutPath"] + "/benchmark/" +
             ProjDirName + "/1_trimming/{porechop_sample}.tsv")
        log:
            (config["OutPath"] + "/logs/" +
             ProjDirName + "/1_trimming/{porechop_sample}.log")
        shell:
            """
            porechop -i {input} \
                     -o {output} \
                     -t {threads}\
                     --extra_end_trim 0 \
                     &> {log}
            """
else:
    if config["Trimmer"] == "cutadapt":  # Choose the trimmer.
        rule all_trim_Illumina:
            input:
                expand(tmp_dir + "/1_trimming/{cutadapt_sample}_R{R_nr}.fq.gz",
                       zip,
                       cutadapt_sample=Samplesheet["sample_name"],
                       R_nr=[1, 2]*len(Samplesheet["sample_name"]))

        rule trim_cutadapt:
            input:
                R1 = config["SamplePath"] + "/" "{cutadapt_sample}_R1.fq.gz",
                R2 = config["SamplePath"] + "/" "{cutadapt_sample}_R2.fq.gz"
            output:
                R1 = tmp_dir + "/1_trimming/{cutadapt_sample}_R1.fq.gz",
                R2 = tmp_dir + "/1_trimming/{cutadapt_sample}_R2.fq.gz"
            threads: config["ThreadNr"]
            conda: "../../envs/preprocessing_env.yaml"
            benchmark:
                (config["OutPath"] + "/benchmark/" +
                 ProjDirName + "/1_trimming/{cutadapt_sample}.tsv")
            log:
                (config["OutPath"] + "/logs/" +
                 ProjDirName + "/1_trimming/{cutadapt_sample}.log")
            shell:
                """
                cutadapt {input} \
                {{'-a ','-A '}}TGAGCTAC{{T{{G{{A,}},}},}}NNNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG\;anywhere\;min_overlap=8\;e=0.1 \
                `#> 5' end adapter; #> brace expansion to generate seq.s; default error tolerance` \
                \
                {{'-g ','-G '}}XAATGATACGGCGACCACCGAGATCTACACNNNNNNNNACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNNNN{{{{{{T,}}C,}}A,}}GTAGCTCA\;anywhere\;min_overlap=8\;e=0.1 \
                `#> 3' end adapter; \
                #> final 8-11 bp are reverse complement of previous seq.; \
                #> X = less strict anchor; default error tolerance` \
                `#> Source: Figure 16 (kit 96D) from SMARTer® ThruPLEX® Tag-seq Kit User Manual, \
                Cat. Nos. R400584, R400585 & R400586 (010719), \
                Takara Bio USA, Inc.` \
                `# Includes "Universal Illumina Adapter".` \
                \
                --times 20 `# ~= read_length / min_overlap` \
                --minimum-length 19 `# (= BWA default minimum seed length)` \
                --cores={threads} \
                -o {output.R1} \
                -p {output.R2} \
                &> {log}
                """

    elif config["Trimmer"] == "bbduk":
        rule all_trimming:
            input:
                expand(tmp_dir + "/1_trimming/{bbduk_sample}_R{R_nr}.fq.gz", zip,
                       bbduk_sample=Samplesheet["sample_name"],
                       R_nr=[1, 2]*len(Samplesheet["sample_name"]))

        rule trim_bbduk:
            input:
                R1 = config["SamplePath"] + "/" "{bbduk_sample}_R1.fq.gz",
                R2 = config["SamplePath"] + "/" "{bbduk_sample}_R2.fq.gz"
            output:
                R1 = tmp_dir + "/1_trimming/{bbduk_sample}_R1.fq.gz",
                R2 = tmp_dir + "/1_trimming/{bbduk_sample}_R2.fq.gz"
            params:
                adapters = config["Adapters"],
                flags = config["BbdukFlags"]
            threads:
                config["ThreadNr"]
            conda: "../../envs/preprocessing_env.yaml"
            benchmark:
                config["OutPath"] + "/benchmark/" +
                ProjDirName + "/1_trimming/{bbduk_sample}.tsv"
            log:
                stats = (config["OutPath"] + "/logs/" +
                         ProjDirName + "/1_trimming/{bbduk_sample}.tsv"),
                Log = (config["OutPath"] + "/logs/" +
                       ProjDirName + "/1_trimming/{bbduk_sample}.log")
            shell:
                """
                bbduk.sh \
                in1={input.R1} \
                in2={input.R2} \
                out1={output.R1} \
                out2={output.R2} \
                ref={params.adapters} \
                {params.flags} \
                stats={log.stats} \
                t={threads} \
                2> {log.Log}
                """

    elif config["Trimmer"] == "none":
        print("\nTrimming is turned off. "
              "To turn it on, please chose a trimmer in the config file!\n")
