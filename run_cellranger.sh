#!/bin/bash
echo "begin"
date
sample=$1
ref=/home/gangcai/shared_space/SoftwareIndex/cellranger/Danio_rerio_gfp/
cellranger count --id "$sample" --sample "$sample" --fastqs="../fq/" --transcriptome $ref --jobmode local --localcores 20 --localmem 70
echo "completed"
date
