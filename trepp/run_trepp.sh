#!/bin/bash

chmod +x feature_tools/trf-4.10.0

bed_file=/home/xuminghua/TREPP/trepp/trepp_test/test.bed
outdir=/home/xuminghua/TREPP/trepp/trepp_test/output
model_dir=/home/xuminghua/TREPP/trepp/models
reference_path=/home/xuminghua/STRP3/code/strp3_2_hg19/data/hs37d5.fa

python predict.py \
 --input_file $bed_file \
 --outdir $outdir \
 --model_dir models \
 --reference $reference_path \
 --db_file genome_annotation.db \
 --suffix trepp \
 --threads 20