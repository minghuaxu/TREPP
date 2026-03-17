#!/bin/bash
python src/trepp_predict.py \
    --ref-fasta demo/hg19_demo.fa \
    --input-bed demo/demo.bed \
    --db-path data/genome_annotations.db \
    --threshold 0.5 \
    --outdir demo/out/predict_run_001