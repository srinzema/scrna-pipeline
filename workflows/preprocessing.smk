from pathlib import Path

# Remove ambient RNA
rule cellbender:
    input: f"{PROJECT_ROOT}/raw_adata/{{sample}}_raw.h5ad"
    output: f"{PROJECT_ROOT}/clean_adata/{{sample}}_denoised.h5"
    log: f"{PROJECT_ROOT}/logs/cellbender/{{sample}}.log"
    params: 
        total_droplets=config["cellbender"]["total_droplets"],
        temp_dir=lambda wildcards: f"/tmp/{wildcards.sample}",
        checkpoint_flag=lambda wildcards: f"--checkpoint /tmp/{wildcards.sample}/ckpt.tar.gz" if Path(f"/tmp/{wildcards.sample}/ckpt.tar.gz").exists() else ""
    resources: 
        nvidia_gpu=1,
        gpu_memory="4G"
    threads: config["cellbender"]["threads"] or 16
    conda: "../envs/cellbender.yaml"
    shell:
        """
        mkdir -p {params.temp_dir};
        cd {params.temp_dir};

        cellbender remove-background --input {input} --output {output} \
        --total-droplets-included 50000 \
        --cpu-threads {threads} {params.checkpoint_flag} \
        --cuda > {log} 2>&1;

        # Remove temp_dir only if it is empty
        rmdir {params.temp_dir} 2>/dev/null || true
        """

# https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html#filtering-low-quality-cells
rule filter_adata:
    input: rules.cellbender.output
    output: f"{PROJECT_ROOT}/filtered_adata/{{sample}}_filtered.h5ad"
    log: f"{PROJECT_ROOT}/logs/filter_adata/{{sample}}.log"
    params: min_genes=config["preprocessing"]["min_genes"] or 200
    threads: 1
    conda: "../envs/filtering.yaml"
    shell: 
        """
        scripts/filter_adata.py \
        --input {input} --output {output} \
        --min_genes {params.min_genes} > {log} 2>&1
        """

# https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html#doublet-detection
rule doublet_detection:
    input: rules.filter_adata.output
    output: f"{PROJECT_ROOT}/adata/{{sample}}_singlet.h5ad"
    log: f"{PROJECT_ROOT}/logs/doublet_detection/{{sample}}.log"
    params:
        n_mads=5,
        n_mads_mt=3
    threads: 1
    conda: "../envs/filtering.yaml"
    shell:
        """
        scripts/doublet_detection.py \
        --input {input} --output {output} \
        --n_mads {params.n_mads} \
        --n_mads_mt {params.n_mads_mt} > {log} 2>&1
        """