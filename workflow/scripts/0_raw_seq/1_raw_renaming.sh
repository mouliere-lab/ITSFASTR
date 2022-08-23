### Rename
for fq in *.fastq.gz; do
    new_filename=$(echo $fq |
                       sed --regexp-extended \
                           --expression='s ([12])\_001.fastq\.gz \1\.fastq\.gz ' \
                           -e 's \.fastq\.gz \.fq\.gz ') \
                           `#> Space as delimiter`
    #echo $new_filename
    mv -v --no-clobber $fq $new_filename
done  # This will change the extention.

for fq in *.fq.gz; do
    new_filename=$(echo $fq |
                       sed --regexp-extended \
                           --expression='s 3\.fq\.gz 2\.fq\.gz ' ) \
                           `#> Space as delimiter`
    #echo $new_filename
    mv -v --no-clobber $fq $new_filename
done
