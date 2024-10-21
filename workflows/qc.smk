from pathlib import Path


rule cellbender:
    input: 
        f"{PROJECT_ROOT}/raw_adata/{{sample}}_raw.h5ad"
    output: 
        f"{PROJECT_ROOT}/clean_adata/{{sample}}_denoised.h5"
    log: 
        f"{PROJECT_ROOT}/logs/cellbender/{{sample}}.log"
    params: 
        temp_dir=lambda wildcards: f"/tmp/{wildcards.sample}",
        checkpoint_flag=lambda wildcards: f"--checkpoint /tmp/{wildcards.sample}/ckpt.tar.gz" if Path(f"/tmp/{wildcards.sample}/ckpt.tar.gz").exists() else ""
    resources: 
        nvidia_gpu=1,
        gpu_memory="4G"
    threads: config["cellbender"]["threads"] or 16
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
    