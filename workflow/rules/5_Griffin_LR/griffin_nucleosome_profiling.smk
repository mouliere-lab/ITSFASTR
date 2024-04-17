#griffin_nucleosome_profiling.snakefile
#Anna-Lisa Doebley
#Template made 2021-12-12
#Ha Lab
#Fred Hutchinson Cancer Research Center

# Modified for ITSFASTR by:
# Norbert Moldovan
# 01-08-2022
# Mouliere Lab
# Amsterdam UMC Cancer Center

import pandas as pd

configfile: "../../../config/config.yaml"  # Set config file.
Samplesheet = pd.read_csv(config["Samplesheet"], sep='\s+')  # Read sample sheet in a dataframe.

# Select output path depending on run mode.
RefPath = config["RefGenome"]
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

rule all:
    input:
        expand("{results_dir}/{samples}/{samples}.GC_corrected.coverage.tsv",
               results_dir=config['OutPath'] + "/" + ProjDirName + "/Griffin_LR",
               samples=Samplesheet["sample_name"]),
        expand("{results_dir}/plots/{site_lists}.pdf",
               results_dir=config['OutPath'] + "/" + ProjDirName + "/Griffin_LR",
               site_lists=config["site_lists"])


rule calc_cov:
    input:
        bam = tmp_dir + "/2_filter_and_index/{samples}.bam",
        GC_bias = config['OutPath'] + "/" + ProjDirName + "/Griffin_LR/GC_bias/{samples}.GC_bias.txt",
    output:
        uncorrected_bw = temp(tmp_dir+"/{samples}/tmp_bigWig/{samples}.uncorrected.bw"),
        GC_corrected_bw = temp(tmp_dir+"/{samples}/tmp_bigWig/{samples}.GC_corrected.bw"),
        tmp_pybedtools = directory(tmp_dir+"/{samples}/tmp_pybedtools")
    params:
        sample_name = "{samples}",
        mappability_bias = 'none',
        mappability_correction = 'False',
        tmp_dir = tmp_dir,
        reference_genome = RefPath,
        mappability_bw = config['mappability_bw'],
        chrom_sizes_path = config['chrom_sizes_path'],
        sites_yaml = "../../../config/config.yaml",
        griffin_scripts_dir = "../../scripts/5_Griffin_LR",
        griffin_coverage_script = "../../scripts/5_Griffin_LR/griffin_coverage.py",
        chrom_column=config['chrom_column'],
        position_column=config['position_column'],
        strand_column=config['strand_column'],
        chroms = config['chroms'],
        norm_window = config['norm_window'],
        size_range=config['size_range'],
        map_quality=config['map_quality'],
        number_of_sites=config['number_of_sites'],
        sort_by=config['sort_by'],
        ascending=config['ascending'],
    threads:  config["ThreadNr"]
    conda: "../../envs/Griffin_LR_env.yaml"
    shell:
        """
		time {params.griffin_coverage_script} \
		--sample_name {params.sample_name} \
		--bam {input.bam} \
		--GC_bias {input.GC_bias} \
		--mappability_bias {params.mappability_bias} \
		--mappability_correction {params.mappability_correction} \
		--tmp_dir {params.tmp_dir} \
		--reference_genome {params.reference_genome} \
		--mappability_bw {params.mappability_bw} \
		--chrom_sizes_path {params.chrom_sizes_path} \
		--sites_yaml {params.sites_yaml} \
		--griffin_scripts {params.griffin_scripts_dir} \
		--chrom_column {params.chrom_column} \
		--position_column {params.position_column} \
		--strand_column {params.strand_column} \
		--chroms {params.chroms} \
		--norm_window {params.norm_window} \
		--size_range {params.size_range} \
		--map_quality {params.map_quality} \
		--number_of_sites {params.number_of_sites} \
		--sort_by {params.sort_by} \
		--ascending {params.ascending} \
		--CPU {threads}
		"""


rule merge_sites:
    input:
        uncorrected_bw = tmp_dir + "/{samples}/tmp_bigWig/{samples}.uncorrected.bw",
        GC_corrected_bw = tmp_dir + "/{samples}/tmp_bigWig/{samples}.GC_corrected.bw",
    output:
        uncorrected_cov = "{results_dir}/{samples}/{samples}.uncorrected.coverage.tsv",
        GC_corrected_cov = "{results_dir}/{samples}/{samples}.GC_corrected.coverage.tsv"
    params:
        sample_name = "{samples}",
        tmp_dir=tmp_dir,
        results_dir="{results_dir}",
        mappability_correction = 'False',
        GC_map_corrected_bw = 'none',
        mappability_bw = config['mappability_bw'],
        chrom_sizes_path = config['chrom_sizes_path'],
        sites_yaml = "../../../config/config.yaml",
        griffin_scripts_dir = "../../scripts/5_Griffin_LR",
        griffin_merge_sites_script = "../../scripts/5_Griffin_LR/griffin_merge_sites.py",
        chrom_column=config['chrom_column'],
        position_column=config['position_column'],
        strand_column=config['strand_column'],
        chroms = config['chroms'],
        norm_window = config['norm_window'],
        save_window = config['save_window'],
        center_window = config['center_window'],
        fft_window = config['fft_window'],
        fft_index = config['fft_index'],
        smoothing_length = config['smoothing_length'],
        encode_exclude = config['encode_exclude'],
        centromeres = config['centromeres'],
        gaps = config['gaps'],
        patches = config['patches'],
        alternative_haplotypes = config['alternative_haplotypes'],
        step = config['step'],
        CNA_normalization = config['CNA_normalization'],
        individual = config['individual'],
        smoothing = config['smoothing'],
        exclude_outliers = config['exclude_outliers'],
        exclude_zero_mappability = config['exclude_zero_mappability'],
        number_of_sites=config['number_of_sites'],
        sort_by=config['sort_by'],
        ascending=config['ascending']
    threads: config["ThreadNr"]
    conda: "../../envs/Griffin_LR_env.yaml"
    shell:
        """
        time {params.griffin_merge_sites_script} \
        --sample_name {params.sample_name} \
        --uncorrected_bw_path {input.uncorrected_bw} \
        --GC_corrected_bw_path {input.GC_corrected_bw} \
        --GC_map_corrected_bw_path {params.GC_map_corrected_bw} \
        --mappability_correction {params.mappability_correction} \
        --tmp_dir {params.tmp_dir} \
        --results_dir {params.results_dir} \
        --mappability_bw {params.mappability_bw} \
        --chrom_sizes_path {params.chrom_sizes_path} \
        --sites_yaml {params.sites_yaml} \
        --griffin_scripts {params.griffin_scripts_dir} \
        --chrom_column {params.chrom_column} \
        --position_column {params.position_column} \
        --strand_column {params.strand_column} \
        --chroms {params.chroms} \
        --norm_window {params.norm_window} \
        --save_window {params.save_window} \
        --center_window {params.center_window} \
        --fft_window {params.fft_window} \
        --fft_index {params.fft_index} \
        --smoothing_length {params.smoothing_length} \
        --exclude_paths {params.encode_exclude} {params.centromeres} {params.gaps} {params.patches} {params.alternative_haplotypes} \
        --step {params.step} \
        --CNA_normalization {params.CNA_normalization} \
        --individual {params.individual} \
        --smoothing {params.smoothing} \
        --exclude_outliers {params.exclude_outliers} \
        --exclude_zero_mappability {params.exclude_zero_mappability} \
        --number_of_sites {params.number_of_sites} \
        --sort_by {params.sort_by} \
        --ascending {params.ascending} \
        --CPU {threads}
        """

rule generate_plots:
    input:
        uncorrected_cov = expand("{results_dir}/{samples}/{samples}.uncorrected.coverage.tsv",
                                 results_dir=config['OutPath'] + "/" + ProjDirName + "/Griffin_LR",
                                 samples=Samplesheet["sample_name"]),
        GC_corrected_cov = expand("{results_dir}/{samples}/{samples}.GC_corrected.coverage.tsv",
                                  results_dir=config['OutPath'] + "/" + ProjDirName + "/Griffin_LR",
                                  samples=Samplesheet["sample_name"]),
    output:
        plots = expand("{results_dir}/plots/{site_lists}.pdf",
                       results_dir=config['OutPath'] + "/" + ProjDirName + "/Griffin_LR",
                       site_lists=config["site_lists"])
    params:
        in_dir = config['OutPath'] + "/" + ProjDirName + "/Griffin_LR",
        samples_yaml = config['OutPath'] + "/" + ProjDirName + "/Griffin_LR/samples.GC.yaml",
        mappability_correction = 'False',
        save_window=config['save_window'],
        step=config['step'],
        individual=config['individual'],
        results_dir=config['OutPath'] + "/" + ProjDirName + "/Griffin_LR",
        griffin_plot_script = "../../scripts/5_Griffin_LR/griffin_plot.py",
        tmp_dir = tmp_dir #for removing empty directories
    conda: "../../envs/Griffin_LR_env.yaml"
    shell:
        """
        time {params.griffin_plot_script} \
		--in_dir {params.in_dir} \
		--samples_yaml {params.samples_yaml} \
		--mappability_correction {params.mappability_correction} \
		--save_window {params.save_window} \
		--step {params.step} \
		--individual {params.individual} \
		--out_dir {params.results_dir}; \
		find {params.tmp_dir} -type d -empty -delete
		"""
