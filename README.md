# TREPP - Tandem Repeat Expansion Pathogenicity Prediction
TREPP is a bioinformatics pipeline for predicting tandem repeat expansions pathogenicity in genomic data.

## Installation

1. Clone this repository:
   ```bash
   git clone https://github.com/minghuaxu/TREPP.git
   cd TREPP/trepp
   ```

2. Ensure dependencies are installed:
    ```bash
   pip install -r requirements.txt
   ```

## Usage
### Basic Command

```bash
    python predict.py \
    --input_file <path_to_bed_file> \
    --outdir <output_directory> \
    --model_dir models \
    --reference <path_to_reference_fasta> \
    --db_file genome_annotation.db \
    --suffix trepp \
    --threads <number_of_threads>
```

## Example
```bash
    python predict.py \
    --input_file /path/to/test.bed \
    --outdir /path/to/output \
    --model_dir models \
    --reference /path/to/hs37d5.fa \
    --db_file genome_annotation.db \
    --suffix trepp \
    --threads 20
```
