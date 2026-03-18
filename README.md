## TREPP: Tandem Repeat Expansion Pathogenicity Prediction

**TREPP** is an interpretable ensemble learning framework designed to robustly identify high-risk TR loci. By combining genomic location information with local sequence properties, TREPP effectively captures the biological context of TR expansions. 

---
### 🛠️ Installation

#### 1. Cloning repository
```bash
git clone https://github.com/minghuaxu/TREPP.git
cd TREPP
```

#### 2. Environment configuration
We recommend using Mamba to manage dependencies for faster environment resolution.
```bash
mamba env create -f trepp_env.yml
conda activate trepp
```
*Note: Please ensure that the `trf` (Tandem Repeats Finder) binary file has execute permissions.*

---
### 📥 Data and Models Download

Before running the prediction process, please manually download the required database files and pre-trained models.

**1. Annotation Database:**
This file must be downloaded for feature extraction.
* 👉 [genome_annotations.db](https://github.com/minghuaxu/TREPP/releases/download/v1.0.0/genome_annotations.db)

**2. Models:**
* 👉 **[trepp_production.pkl](https://github.com/minghuaxu/TREPP/releases/download/v1.0.0/trepp_production.pkl)**: This is a production environment model used for real-world forecasting. It was trained using **all available data** and exhibits optimal generalization ability and predictive performance.
* 👉 [trepp_experiments.pkl](https://github.com/minghuaxu/TREPP/releases/download/v1.0.0/trepp_experiments.pkl): This is the experimental environment model. It was trained using only the training set data and is primarily used to reproduce the cross-validation and independent test set evaluation results presented in our research paper.

**Example of a quick download command:**
```bash
# Download the database file to the data directory
wget https://github.com/minghuaxu/TREPP/releases/download/v1.0.0/genome_annotations.db -P data/

# Download the production model to the models directory
wget https://github.com/minghuaxu/TREPP/releases/download/v1.0.0/trepp_production.pkl -P src/models/
```

---
### 📄 Input Format

TREPP accepts standard headerless, tab-separated BED files as input.

The first four columns of the file must be arranged strictly in the following order:
1. **Chromosome**: Chromosome number (e.g., `chr1`, `chrX`)
2. **Start**: Start coordinates of the TR interval (0-based)
3. **End**: End coordinates of the TR interval (1-based)
4. **Motif**: Base sequence of the tandem repeat unit (e.g., `CAG`, `ATCC`)

**Input Example (`example.bed`):**
```text
chr4	3076604	3076667	CAG
chr12	306240	306372	ATCC
```
*💡 Tip: If your BED file contains 5 or more columns (e.g., actual labels or other annotations), the program will automatically ignore the extra columns and extract only the first four for calculation.*

---
### 📖 Usage

#### Help
```bash
python trepp_predict.py -h
```

#### Predict directly from BED
```bash
python src/trepp_predict.py \
  --ref-fasta example/hg19_example.fa \
  --input-bed example/example.bed \
  --db-path data/genome_annotations.db \
  --outdir example/out/predict_run_001
```

Artifacts produced:
- `out/predict_run_001/intermediate/predict_features.csv`
- `out/predict_run_001/logs/extract_features.log`
- `out/predict_run_001/trepp_predicted.tsv`

#### Override prediction output path
If you want a custom filename, still keep it inside the run folder:
```bash
python src/trepp_predict.py \
  --ref-fasta example/hg19_example.fa \
  --input-bed example/example.bed \
  --db-path data/genome_annotations.db \
  --outdir out/predict_run_002 \
  --pred-out out/predict_run_002/custom_predictions.tsv
```

#### Thresholding
If `--threshold` is omitted, the script defaults to use `0.5`. 
Custom threshold:
```bash
python trepp_predict.py \
  --ref-fasta example/hg19_example.fa \
  --input-bed example/example.bed \
  --db-path data/genome_annotations.db \
  --outdir out/predict_run_003 \
  --threshold 0.35
```

---
## 📊 Output format

The prediction output is a TSV containing:
- `id` (from input `id`, or constructed from `chr/start/end/motif`, else row index)
- `prob` (final stacked probability)
- `pred_label` (`prob >= threshold`)

---
## 📜 Citation
If you used TREPP in your research, please cite it:
```bibtex
Xu M, Hu K, Deng C, et al. TREPP: Tandem Repeat Expansion Pathogenicity Predicting Approach Using Stacked CatBoost Models and Multiple Features[C]//International Symposium on Bioinformatics Research and Applications. Singapore: Springer Nature Singapore, 2025: 325-336.
```