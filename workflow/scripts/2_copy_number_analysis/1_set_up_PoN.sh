# this script describes how to set up panel of normals
## NOTE: You don't need to generate a new panel of normals, there are two reference sets in /projects/0/fmlab/ref_genomes/ref_ichorCNA_PoN/
## In the ichorCNA.smk, the panel of normals generated from STM SLX-13222 is automatically included

# after generating the .wig files using readcounter script (implemented in snakemake):

# generate PoN by running the createPanelOfNormals.R script
Rscript ichorCNA/scripts/createPanelOfNormals.R \
    --filelist healthy_ref_wig_paths.txt \
    --gcWig ichorCNA/inst/extdata/gc_hg38_1000kb.wig --mapWig ichorCNA/inst/extdata/map_hg38_1000kb.wig \
    --centromere ichorCNA/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt \
    --outfile panel_of_normals

# --filelist /path/to/wig_files.txt - file containing a list of the paths to all the normals in the panel to analyze; each listed on a separate line
