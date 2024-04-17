# griffin_GC_and_mappability_correction.smk
# Original code by:
# Anna-Lisa Doebley
# 2021-12-13
# Ha Lab
# Fred Hutchinson Cancer Research Center

# Modified for ITSFASTR by:
# Norbert Moldovan
# 01-08-2022
# Mouliere Lab
# Amsterdam UMC Cancer Center

import pandas as pd

configfile: "../../../config/config.yaml"  # Set config file.
Samplesheet = pd.read_csv(config["Samplesheet"], sep='\s+')  # Read sample sheet in a dataframe.

# Select output path depending on run mode.
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

if config['mappability_correction']:
    rule all:
        input:
            expand("{out_dir}/mappability_bias/{samples}.mappability_bias.txt",
                samples=Samplesheet["sample_name"],
                out_dir=config['OutPath'] + "/" + ProjDirName + "/Griffin_LR"),
            expand("{out_dir}/mappability_plots/{samples}.mappability_bias.pdf",
                samples=Samplesheet["sample_name"],
                out_dir=config['OutPath'] + "/" + ProjDirName + "/Griffin_LR"),
            expand("{out_dir}/mappability_plots/{samples}.mappability_bias.read_coverage_distribution.pdf",
                samples=Samplesheet["sample_name"],
                out_dir=config['OutPath'] + "/" + ProjDirName + "/Griffin_LR"),
            expand("{out_dir}/GC_counts/{samples}.GC_counts.txt",
                samples=Samplesheet["sample_name"],
                out_dir=config['OutPath'] + "/" + ProjDirName + "/Griffin_LR"),
            expand("{out_dir}/GC_bias/{samples}.GC_bias.txt",
                samples=Samplesheet["sample_name"],
                out_dir=config['OutPath'] + "/" + ProjDirName + "/Griffin_LR"),
            expand("{out_dir}/GC_plots/{samples}.GC_bias.summary.pdf",
                samples=Samplesheet["sample_name"],
                out_dir=config['OutPath'] + "/" + ProjDirName + "/Griffin_LR"),
            expand("{out_dir}/samples.GC.yaml",
                out_dir=config['OutPath'] + "/" + ProjDirName + "/Griffin_LR")

else:
    rule all:
        input:
            expand("{out_dir}/GC_counts/{samples}.GC_counts.txt",
                samples=Samplesheet["sample_name"],
                out_dir=config['OutPath'] + "/" + ProjDirName + "/Griffin_LR"),
            expand("{out_dir}/GC_bias/{samples}.GC_bias.txt",
                samples=Samplesheet["sample_name"],
                out_dir=config['OutPath'] + "/" + ProjDirName + "/Griffin_LR"),
            expand("{out_dir}/GC_plots/{samples}.GC_bias.summary.pdf",
                samples=Samplesheet["sample_name"],
                out_dir=config['OutPath'] + "/" + ProjDirName + "/Griffin_LR"),
            expand("{out_dir}/samples.GC.yaml",
                out_dir=config['OutPath'] + "/" + ProjDirName + "/Griffin_LR")


rule filter_and_index_Griffin_LR:
    input:
        config["OutPath"] + "/" + ProjDirName + "/1_mapping/{samples}.bam"
    output:
        bam = tmp_dir + "/2_filter_and_index/{samples}.bam",
        bai = tmp_dir + "/2_filter_and_index/{samples}.bam.bai"
    threads: config["ThreadNr"]
    params:
        qual = config["FiltQual"]
    conda: "../../envs/Griffin_LR_env.yaml"
    benchmark:
        (config["OutPath"] + "/benchmark/" +
         ProjDirName + "/2_filter_and_index/{samples}.tsv")
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

if config['mappability_correction']:
    # If mappability correction is applied.
    rule mappability_bias:
        input:
            bam_file = tmp_dir + "/2_filter_and_index/{samples}.bam"
        output:
            mappability_bias = "{out_dir}/mappability_bias/{samples}.mappability_bias.txt",
            mappability_plot = "{out_dir}/mappability_plots/{samples}.mappability_bias.pdf",
            mappability_plot2 = "{out_dir}/mappability_plots/{samples}.mappability_bias.read_coverage_distribution.pdf",
            tmp = temp(directory(tmp_dir + "/Griffin_LR/tmp_{samples}/"))
        threads: config["ThreadNr"]
        params:
            sample_name = "{samples}",
            mappability_bw = config['mappability_bw'],
            encode_exclude = config['encode_exclude'],
            centromeres = config['centromeres'],
            gaps = config['gaps'],
            patches = config['patches'],
            alternative_haplotypes = config['alternative_haplotypes'],
            chrom_sizes = config['chrom_sizes'],
            map_q=config['FiltQual'],
            griffin_mappability_script = "../../scripts/5_Griffin_LR/griffin_mappability_correction.py"
        conda: "../../envs/Griffin_LR_env.yaml"
        shell:
            """
            time {params.griffin_mappability_script} \
            --bam_file {input.bam_file} \
            --bam_file_name {params.sample_name} \
            --output {output.mappability_bias} \
            --output_plot {output.mappability_plot} \
            --mappability {params.mappability_bw} \
            --exclude_paths {params.encode_exclude} {params.centromeres} {params.gaps} {params.patches} {params.alternative_haplotypes} \
            --chrom_sizes {params.chrom_sizes} \
            --map_quality {params.map_q} \
            --CPU {threads} \
            --tmp_dir {output.tmp}
            """

rule GC_counts:
    input:
        bam_file = tmp_dir + "/2_filter_and_index/{samples}.bam"
    output:
        GC_counts_file = protected("{out_dir}/GC_counts/{samples}.GC_counts.txt"),
    params:
        sample_name = "{samples}",
        output_dir = "{out_dir}",
        mappable_regions_path = config['mappable_regions'],
        ref_seq = config['RefGenome'],
        chrom_sizes = config['chrom_sizes'],
        map_q=config['FiltQual'],
        size_range=config['GC_bias_size_range'],
        griffin_GC_counts_script = "../../scripts/5_Griffin_LR/griffin_GC_counts.py"
    threads: config["ThreadNr"]
    conda: "../../envs/Griffin_LR_env.yaml"
    shell:
        """
        time {params.griffin_GC_counts_script} \
        --bam_file {input.bam_file} \
        --bam_file_name {params.sample_name} \
        --mappable_regions_path {params.mappable_regions_path} \
        --ref_seq {params.ref_seq} \
        --chrom_sizes {params.chrom_sizes} \
        --out_dir {params.output_dir} \
        --map_q {params.map_q} \
        --size_range {params.size_range} \
        --CPU {threads}
        """

rule GC_bias:
    input:
        GC_counts_file = "{out_dir}/GC_counts/{samples}.GC_counts.txt"
    output:
        GC_bias_file = "{out_dir}/GC_bias/{samples}.GC_bias.txt",
        GC_plots_file = "{out_dir}/GC_plots/{samples}.GC_bias.summary.pdf"
    params:
        sample_name = "{samples}",
        output_dir = "{out_dir}",
        mappable_name =  config['mappable_regions'].rsplit('/',1)[1].rsplit('.',1)[0],
        genome_GC_frequency = config['genome_GC_frequency'],
        size_range=config['GC_bias_size_range'],
        griffin_GC_bias_script =  "../../scripts/5_Griffin_LR/griffin_GC_bias.py"
    conda: "../../envs/Griffin_LR_env.yaml"
    shell:
        """
        time {params.griffin_GC_bias_script} \
        --bam_file_name {params.sample_name} \
        --mappable_name {params.mappable_name} \
        --genome_GC_frequency {params.genome_GC_frequency} \
        --out_dir {params.output_dir} \
        --size_range {params.size_range}
        """

if config['mappability_correction']:
    rule make_samples_yaml:
        input:
            mappability_bias = expand("{out_dir}/mappability_bias/{samples}.mappability_bias.txt",
                samples=Samplesheet["sample_name"],
                out_dir=config['OutPath'] + "/" + ProjDirName + "/Griffin_LR"),
            GC_bias_file = expand("{out_dir}/GC_bias/{samples}.GC_bias.txt") # Needs to exist for the abspath to work.
        output:
            samples_yaml = "{out_dir}/samples.GC.yaml"
        params:
            sample_dict = Samplesheet["sample_name"],
        run:
            import os
            with open(config['OutPath'] + "/" + ProjDirName + "/Griffin_LR/samples.GC.yaml", 'w') as f:
                f.write('samples:\n')
                for key in params['sample_dict']:
                    f.write('  '+key+':\n')
                    f.write('    bam: ' + tmp_dir + "/2_filter_and_index/{samples}.bam\n")
                    GC_bias_path = config['OutPath']+'/Griffin_LR/GC_bias/'+key+'.GC_bias.txt'
                    GC_bias_path = os.path.abspath(GC_bias_path)
                    f.write('    GC_bias: '+GC_bias_path+'\n')
                    mappability_bias_path = config['OutPath'] + ProjDirName + '/Griffin_LR/mappability_bias/'+key+'.mappability_bias.txt'
                    mappability_bias_path = os.path.abspath(mappability_bias_path)
                    f.write('    mappability_bias: '+mappability_bias_path+'\n')
else:
    rule make_samples_yaml:
        input:
            GC_bias_file = expand("{out_dir}/GC_bias/{samples}.GC_bias.txt",
                                  samples=Samplesheet["sample_name"],
                                  out_dir=config['OutPath'] + "/" + ProjDirName + "/Griffin_LR")  #Needs to exist for the abspath to work.
        output:
            samples_yaml = "{out_dir}/samples.GC.yaml"
        params:
            sample_list = Samplesheet["sample_name"],
        run:
            import os
            with open(config['OutPath'] + "/" + ProjDirName + "/Griffin_LR/samples.GC.yaml", 'w') as f:
                f.write('samples:\n')
                for key in Samplesheet["sample_name"]:
                    f.write('  ' + key + ':\n')
                    f.write('    bam: ' + tmp_dir + "/2_filter_and_index/" + key + ".bam\n")
                    GC_bias_path = config['OutPath'] + "/" + ProjDirName + "/Griffin_LR/GC_bias/" + key + ".GC_bias.txt"
                    GC_bias_path = os.path.abspath(GC_bias_path)
                    f.write('    GC_bias: '+GC_bias_path+'\n')