# ITSFASTR - InTegrated Sequence and Fragmentome AnalysiS Time Reduction

## Contents
- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Demo](#demo)
- [Results](#results)
- [License](./LICENSE)
- [Citation](#citation)

## Overview
During cell death the DNA of tumor cells is released in to the bloodstream,
forming a pool of circulating cell-free DNA (cfDNA). cfDNA can be isolated 
from blood and/or urine of cancer patients and analyzed with sequencing. 
Most conventional short-read sequencing methods are technically challenging,
labor intensive and time consuming, requiring several days but more typically 
weeks to obtain interpretable data which are limited by a bias for short cfDNA fragments.
The ITSFASTR toolkit provides the bioinformatic workflow for the analysis of copy number
aberrations, fragmentation and fragment end sequences and nucleosome positioning from 
conventional Illumina and Oxford Nanopore Technologies (ONT) cfDNA data.

## Repo Contents
 - [ITSFASTR pipeline](./workflow)

## System Requirements

### Hardware Requirements
The ITSFASTR pipeline is designed to run on HPC clusters, as some steps can be
resource-heavy with real-world data sets.

If ran on a standard computer, we recommend:
RAM: 16+ GB
CPU: 4+ cores, 3.3+ GHz/core

### Software Requirements

The package is tested on an Ubuntu 20.04.4 LTS operating system, but it should
work on [Windows Subsystem for Linux (WSL)][1] as well.

## Installation Guide

### Package Dependencies
Prior to downloading ITSFASTR users should install [Anaconda][2].

ITSFASTR uses [Snakemake][3] for reproducibility.
Install Snakemake by running:
```
conda create -n snakemake -c bioconda snakemake
```
Activate your environment:
```
conda activate snakemake
```
Snakemake uses Mamba, a fast package manager. To install it type:
```
conda install mamba -n base -c conda-forge
```
Snakemake requires strict channel priorities in Conda, which is also good for
reproducibility. To set this, type:
```
conda config --set channel_priority strict
```
Every other dependency will be installed automatically during the first run of
the pipeline.
Software versions used by the pipeline can be found in the environment files
in the [envs](./workflow/envs/) folder.

In short the pipeline uses:
[BBduk][6] or [Cutadapt][7] for adapter trimming of Illumina reads or 
[Porechop][14] for ONT reads.
[BWA MEM][8] for mapping Illumina and [Minimap2][15] for ONT reads.
[Samtools][9] and [Sambamba][10] for quality filtering and marking duplicates.
[Qualimap][11] for quality control summarized by [MultiQC][12].

### Downloading ITSFASTR
Download ITSFASTR by:
```
wget https://github.com/mouliere-lab/ITSFASTR/archive/refs/heads/main.zip &&
unzip main.zip &&
mv ITSFASTR-main ITSFASTR &&
rm main.zip
```
or if you have Git installed:
```
git clone https://github.com/mouliere-lab/ITSFASTR.git
```

## Demo
The ITSFASTR pipeline accepts paired Illumina sequencing files with 
`.fq.gz` extension. If your files have a different extension use the 
script in `workflow/scripts/0_raw_seq/1_raw_renaming.sh` to change it.\
ONT sequencing files can be fed into the pipeline after demultiplexing.
For an ONT run a folder containing the sequencing files is the input,
for ex. each barcode folder. For further explanation see the structure of
the provided dummy data.

### Dummy Data Set
To run this demo please download [these dummy files][4].
The dummy files are a synthetic admixtures of real Illumina and ONT 
sequencing files from the plasma of cancer patients, downsampled to 100K reads.

### Creating the Sample Sheet
The [sample sheet](./config/samplesheet.csv) is a space-delimited `.csv` file
with two labels in its header: `sample_name group`

An example containing the dummy file names can be found in the
[config](./config) directory.

`sample_name` is the unique name of the sample (file). Paired input files 
have one unique sample name, which is the file name minus
the `_R1(2).fq.gz` suffix.

`group` represents the single analysis group the given file belongs to. For the
final steps of the FrEIA analysis a control group is needed. Mark this with a
`*` following the group name.

### Using the Config File
The [config file](./config/config.yaml) is a `.yaml` format file containing the
parameters used during the run.

#### 1. Set the project name
This will be used as the name of the main folder containing results.
#### 2. Set the run type
Based on the sequencing files either select Illumina or ONT.
#### 3. Set the sample path
For an Illumina run add the path to the folder containing the `.fq.gz` files.\
For an ONT run add the path to the folder containing the output folder(s) of
base calling or demultiplexing.
#### 4. Set the path to the sample sheet
Add the path to the sample sheet file.
#### 5. Set the path to the indexed reference genome
The pipeline maps the reads to a reference genome using either BWA MEM for
Illumina runs or Minimap2 for ONT runs. To map the reads you first need
to index the reference using BWA MEM ([tutorial][8]) for Illumina runs or Minimap2 ([tutorial][16]) for ONT
runs. Indexing large genomes can take from minutes to hours
butt needs to be done only once.
Add the path to the indexed reference genome based on the run type.
#### 6. Set the output path
Add the path to the output folder. A folder with the project name will be
created here containing all the results.
#### 7. Set the path to a temporary folder
The temporary folder (tmp) will be used to store temporary files. On HPC clusters
this should be on a `scratch drive` with fast read/write capability.
#### 8. Set the number of CPU cores the pipeline can use
The number of CPU cores should be less or equal to the available cores.
On a normal computer this is usually 4 or 8.
Be aware that with the increase of the core number the memory usage can 
also increase.
#### 9. Download the mappability track
The mappability track used by Griffin LR is big, thus it is not provided.
[Download][17] a mappability track and place it in [/resources/Griffin_LR](resources/Griffin_LR).
#### 10. Set path to the reference genome
Set the path to the reference genome fasta file.
#### 11. Create the sites of interest file
Griffin LR requires one or more sites of interest file. Each file can contain multiple sites of the same tipe. For an example please see the [Running ITSFASTR on a local machine using Snakemake](#running-itsfastr-on-a-local-machine-using-snakemake).

### Running ITSFASTR on a HPC cluster using Snakemake and Slurm
This is the recommended way of running the pipeline, as some steps are
resource-heavy.
Rules of the pipeline are run separately. To run a rule, first change directories
into the rule's folder for example `cd workflow/rules/1_preprocessing/`
and run snakemake with the `-s` argument followed by the rule name

For example trim the reads by running:
```
snakemake --printshellcmds \
          --keep-going \
          --use-conda \
          --cluster "sbatch \
                     --time=60 \
                     --nodes=1 \
                     --partition=thin \
                     --cpus-per-task=32 \
                     --mem=64000 \
                     --mail-type=FAIL,TIME_LIMIT \
                     --mail-user=your@email.address \
                     --output=path/to/slurm_out/slurm-%j.out" \
          --jobs 500 \
          --max-jobs-per-second 3 \
          --max-status-checks-per-second 5 \
          -s trimming.smk
```
### Running ITSFASTR on a local machine using Snakemake
Running on a local machine is not recommended as some steps may exceed available
resources.
Rules of the pipeline are run separately. To run a rule, first change
directories into the rule's folder and run snakemake with the `-s` argument
followed by the rule name.

In the Demo we are running ITSFASTR on a local machine.
Expected run time: ~ 1.5 h

First download the [human reference genome][18] (GRCh38_no_alt_analysis_set) and
the [mappability track][17]. Extract the reference genome, and move the mappability
track to [/resources/Griffin_LR](resources/Griffin_LR).

For the nucleosome positioning analysis we use nucleosome covered quiescent genomic sites (Quies)
and nuclosome free transcriptional start sites (TssA). Download files containg the genomic positions
of these sites from [here][19].
Set the path to the files containing the sites of interest in the config file `site_list`.

#### 1.1. Trimming
Change directory: `cd workflow/rules/1_preprocessing/`

Remove sequencing adapters by running:
```
snakemake --printshellcmds \
          --keep-going \
          --use-conda \
          --cores 8 \
          -s trimming.smk
```

#### 1.2. Mapping (non-xenograft samples)
Align the reads to the previously set reference genome by running:
```
snakemake --printshellcmds \
          --keep-going \
          --use-conda \
          --cores 8 \
          -s mapping.smk
```


#### 2. Run ichorCNA copy number analysis
Change directory: `cd ../2_copy_number_analysis`.  
The ichorCNA tool requires the path to the ichorCNA script to be set in the config file. To set the path first run:

```
snakemake --printshellcmds \
          --keep-going \
          --use-conda \
          --cores 8 \
          -s ichorCNA.smk \
          --conda-create-envs-only
```
This will create the environment and provide you with a path to it, similar to `.snakemake/conda/103fe64931ac3076f3305793328dd331_`. Paste this path into the [config file](./config/config.yaml) after the `conda_dir:` parameter.

Now run:
```
snakemake --printshellcmds \
          --keep-going \
          --use-conda \
          --cores 8 \
          -s ichorCNA.smk
```
#### 3. Run the fragment size analysis
Change directory: `cd ../3_fragmentation`.

```
snakemake --printshellcmds \
          --keep-going \
          --use-conda \
          --cores 8 \
          -s fragmentation.smk
```

#### 4.1. Run FrEIA fragment end sequence analysis
Change directory: `cd ../4_FrEIA`.

Extract the fragment ends by running:
```
snakemake --printshellcmds \
          --keep-going \
          --use-conda \
          --cores 8 \
          -s preprocessing.smk
```

#### 4.2. Calculate fragment end proportions and diversity with FrEIA:
```
snakemake --printshellcmds \
          --keep-going \
          --use-conda \
          --cores 8 \
          -s FrEIA.smk
```

#### 5.1. Compute GC content correction with Griffin LR:
Change directory: `cd ../5_Griffin_LR`.

```
snakemake --printshellcmds \
          --keep-going \
          --use-conda \
          --cores 8 \
          -s griffin_GC_and_mappability_correction.smk
```

#### 5.2. Perform nucleosome profiling using Griffin LR:

```
snakemake --printshellcmds \
          --keep-going \
          --use-conda \
          --cores 8 \
          -s griffin_nucleosome_profiling.smk
```

### Results
Please download the expected results from [here][13]. The expected results
archive contains the expected folder structure, quality check results, the copy number
profile, the fragment length profile, the final files fragment end sequence proportions, diversity 
and the nucleosome profiles of the dummy data set.

### Xenomapping (xenograft samples only)
If starting from xenograft samples reads from the host and the graft need to be separated after trimming. To do this, first [create the index files for both the host and the graft genomes](#5.-set-the-path-to-the-indexed-reference-genome). \
Then add the paths to the [config file](./config/config.yaml) after `RefPath_Primary` and `RefPath_Secondary`. \
Finally change directory: `cd workflow/rules/6_xenomapping/`
and run:
```
snakemake --printshellcmds \
          --keep-going \
          --use-conda \
          --cores 8 \
          -s xenomapping_LR.smk
```
The deconvoluted and mapped reads are stored in the same folder as the output of the normal mapping is.

## Citation
For usage of the ITSFASTR pipeline please cite:
van der Pol et al. Real-time analysis of the cancer genome and fragmentome from plasma and urine short and long cell-free DNA using Nanopore sequencing (2022)

[1]: https://docs.microsoft.com/en-us/windows/wsl/install
[2]: https://docs.anaconda.com/anaconda/install/linux/
[3]: https://snakemake.github.io/
[4]: https://doi.org/10.6084/m9.figshare.20550741
[5]: https://github.com/epigenelabs/pyComBat
[6]: https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/
[7]: https://cutadapt.readthedocs.io/en/stable/
[8]: http://bio-bwa.sourceforge.net/bwa.shtml
[9]: http://www.htslib.org/
[10]: https://lomereiter.github.io/sambamba/
[11]: http://qualimap.conesalab.org/
[12]: https://multiqc.info/
[13]: https://doi.org/10.6084/m9.figshare.20550960.v1
[14]: https://github.com/rrwick/Porechop
[15]: https://github.com/lh3/minimap2
[16]: https://lh3.github.io/minimap2/minimap2.html
[17]: https://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/k50.Umap.MultiTrackMappability.bw
[18]: http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/analysisSet/hg38.fullAnalysisSet.chroms.tar.gz
[19]: https://doi.org/10.6084/m9.figshare.20551620.v1