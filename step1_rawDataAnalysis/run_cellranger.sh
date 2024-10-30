#!/bin/bash
echo "begin"
date
sample=${arg1}
echo "start count for "$sample
date
ref=/home/db/public/annotation/SoftwareIndex/2019Version/cellranger/Danio_rerio/Danio_rerio_gfp/
fq="../combined/"
cellranger count --jobmode local --localcores 15 --localmem 60 --id "$sample" --sample "$sample" --fastqs $fq --transcriptome $ref
echo "end count for $sample"
date
