configfile: "config/config.yaml"

rule all:
    input:
        "results/performance_benchmark.csv",
        "results/loco_performance.csv",
        "results/ablation_study.csv" 

# ---------------------------------------------------------------------
# 1. 构建数据库 (Build Database)
# ---------------------------------------------------------------------
rule build_database:
    input:
        gff = config["annotations"]["gff"],
    output:
        db = config["processed"]["db"]
    params:
        script = "scripts/build_database.py"
    shell:
        "python {params.script} --gff {input.gff} --gnomad {input.gnomad} --out {output.db}"

# ---------------------------------------------------------------------
# 2. 特征提取 (Unified Feature Extraction)
# ---------------------------------------------------------------------
# 现在只需要输入 BED, FASTA 和 DB
rule extract_features_train:
    input:
        bed = config["datasets"]["train"],
        ref = config["reference"]["fasta"],
        db = config["processed"]["db"],
        trf = config["paths"]["trf_binary"]  # 确保这是一个文件路径
    output:
        features = config["processed"]["train_features"]
    params:
        script = "scripts/extract_features.py"
    log:
        "logs/extract_features_train.log"
    shell:
        """
        # 给 TRF 加执行权限，防止 Permission denied
        chmod +x {input.trf}
        
        python {params.script} \
            --input_bed {input.bed} \
            --ref_fasta {input.ref} \
            --db_path {input.db} \
            --trf_path {input.trf} \
            --output_csv {output.features} \
            > {log} 2>&1
        """

rule extract_features_eval:
    input:
        bed = config["datasets"]["eval"],
        ref = config["reference"]["fasta"],
        db = config["processed"]["db"],
        trf = config["paths"]["trf_binary"]
    output:
        features = config["processed"]["eval_features"]
    params:
        script = "scripts/extract_features.py"
    log:
        "logs/extract_features_eval.log"
    shell:
        """
        # 给 TRF 加执行权限，防止 Permission denied
        chmod +x {input.trf}
        
        python {params.script} \
            --input_bed {input.bed} \
            --ref_fasta {input.ref} \
            --db_path {input.db} \
            --trf_path {input.trf} \
            --output_csv {output.features} \
            > {log} 2>&1
        """

rule analyze_features:
    input:
        features = config["processed"]["train_features"]
    output:
        plot1 = "results/plots/plot_clustering.png",
        plot2 = "results/plots/plot_importance.png"
    params:
        script = "scripts/analyze_separability.py",
        out_dir = "results/plots"
    shell:
        """
        python {params.script} \
            --input {input.features} \
            --out_dir {params.out_dir}
        """

# ---------------------------------------------------------------------
# 3. 模型训练 (Model Training & New Experiments)
# ---------------------------------------------------------------------
rule run_experiments:
    input:
        train_csv = config["processed"]["train_features"],
        eval_csv = config["processed"]["eval_features"]
    output:
        bench = "results/performance_benchmark.csv",
        loco = "results/loco_performance.csv",
        probs_npy = "results/probs_TREPP_(Ours).npy"
    params:
        script = "scripts/train_advanced.py",
        out_dir = "results",
        model_dir = "data/models",
        n_models = 15,  # Number of base models
        n_trials = 150  # Optuna trials (increased for better tuning)
    log:
        "logs/train_advanced.log"
    shell:
        """
        python {params.script} \
            --train {input.train_csv} \
            --eval {input.eval_csv} \
            --out_dir {params.out_dir} \
            --model_dir {params.model_dir} \
            --n_models {params.n_models} \
            --n_trials {params.n_trials} \
            > {log} 2>&1
        """

# 定义要遍历的模型数量范围 (7 到 14)
N_MODELS_RANGE = range(2, 11)

# 1. 汇总规则：请求所有结果文件
rule run_model_sweep:
    input:
        expand("results/performance_benchmark_{n}.csv", n=N_MODELS_RANGE),

# 2. 动态运行规则
rule run_experiments_by_n:
    input:
        train_csv = config["processed"]["train_features"], 
        eval_csv = config["processed"]["eval_features"]
    output:
        bench = "results/performance_benchmark_{n}.csv",
        model = "data/models_{n}/trepp_final_{n}.pkl"
    params:
        script = "scripts/trepp_v7_macro.py",
        # 创建唯一的临时目录，防止并行冲突
        temp_out_dir = "results/temp_run_{n}",
        model_dir = "data/models_{n}",
        n_models = "{n}",
        n_trials = 150 # 扫参时可以稍微减少trial次数，或者保持150
    log:
        "logs/train_advanced_{n}.log"
    threads: 24
    shell:
        """
        # 1. 创建临时输出目录
        mkdir -p {params.temp_out_dir}

        # 2. 运行 Python 脚本
        # 注意：这里我们让脚本输出到临时目录
        python {params.script} \
            --train {input.train_csv} \
            --eval {input.eval_csv} \
            --out_dir {params.temp_out_dir} \
            --model_dir {params.model_dir} \
            --n_models {params.n_models} \
            --n_trials {params.n_trials} \
            > {log} 2>&1

        # 3. 将结果文件重命名并移至最终位置
        mv {params.temp_out_dir}/performance_benchmark.csv {output.bench}
        
        # 4. 重命名并移动模型文件 (假设脚本默认保存为 trepp_final.pkl)
        mv {params.model_dir}/trepp_final.pkl {output.model}

        # 5. 清理临时目录
        rm -rf {params.temp_out_dir}
        """


# 新增消融实验规则
rule run_ablation:
    input:
        train_csv = "data/processed/train_features.csv",
        eval_csv = "data/processed/eval_features.csv"
    output:
        csv = "results/ablation_study.csv"
    params:
        script = "scripts/run_ablation.py"
    log:
        "logs/ablation.log"
    shell:
        """
        python {params.script} \
            --train {input.train_csv} \
            --eval {input.eval_csv} \
            --out_csv {output.csv} \
            > {log} 2>&1
        """

