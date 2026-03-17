#!/bin/bash
python src/trepp_predict.py \
    --ref-fasta example/hg19_example.fa \
    --input-bed example/example.bed \
    --db-path data/genome_annotations.db \
    --threshold 0.5 \
    --outdir example/out/predict_run_001