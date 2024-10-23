configfile: "config.yaml"
PROJECT_ROOT = config["project_root"]

include: "workflows/preprocessing.smk"
micromamba = True

samples = [x.name.strip("_raw.h5ad") for x in Path(f"{PROJECT_ROOT}/raw_adata").iterdir()]

rule all:
    input: 
        expand(f"{PROJECT_ROOT}/clean_adata/{{sample}}_denoised.h5", sample=samples), # cellbender
        expand(f"{PROJECT_ROOT}/filtered_adata/{{sample}}_filtered.h5", sample=samples)