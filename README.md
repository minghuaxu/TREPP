## TREPP: Tandem Repeat Evaluation & Pathogenicity Predictor
![alt text](https://img.shields.io/badge/workflow-Snakemake-brightgreen.svg)
![alt text](https://img.shields.io/badge/License-MIT-yellow.svg)
![alt text](https://img.shields.io/badge/python-3.9+-blue.svg)


This package provides an **end-to-end prediction** entrypoint for TREPP:
1) **feature extraction** from BED via your existing `scripts/extract_features.py`
2) **feature alignment** to the training feature list stored in `trepp_final.pkl`
3) **stacked inference** (Calibrated base models → Logistic Regression meta model)

All intermediate files (features, logs) and final outputs are written under a user-specified `--outdir`.

### 🛠️ 安装指南
---
#### 1. 克隆仓库
```Bash
git clone https://github.com/minghuaxu/TREPP.git
cd TREPP
```

#### 2. 环境配置
推荐使用 Conda 或 Mamba 管理依赖：
```Bash
conda env create -f environment.yaml
conda activate trepp
```
注意：请确保 trf (Tandem Repeats Finder) 二进制文件已获得执行权限。

---
### 📖 Usage

#### Help

```bash
python trepp_predict.py -h
```

### Predict directly from BED

```bash
python trepp_predict.py \
  --input-bed data/predict.bed \
  --ref-fasta data/reference.fa \
  --outdir out/predict_run_001
```

Artifacts produced:
- `out/predict_run_001/intermediate/predict_features.csv`
- `out/predict_run_001/logs/extract_features.log`
- `out/predict_run_001/predictions.tsv`

### Override prediction output path

If you want a custom filename, still keep it inside the run folder:

```bash
python trepp_predict.py \
  --input-bed data/predict.bed \
  --ref-fasta data/reference.fa \
  --outdir out/predict_run_002 \
  --pred-out out/predict_run_002/custom_predictions.tsv
```

### Thresholding

If `--threshold` is omitted, the script default use `0.5`.

Custom threshold:

```bash
python trepp_predict.py \
  --input-bed data/predict.bed \
  --ref-fasta data/reference.fa \
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
## 📜 引用
如果您在研究中使用了 TREPP，请引用：
```bash
Xu M, Hu K, Deng C, et al. TREPP: Tandem Repeat Expansion Pathogenicity Predicting Approach Using Stacked CatBoost Models and Multiple Features[C]//International Symposium on Bioinformatics Research and Applications. Singapore: Springer Nature Singapore, 2025: 325-336.
```
