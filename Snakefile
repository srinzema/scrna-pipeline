configfile: "config.yaml"
PROJECT_ROOT = config["project_root"]

include: "workflows/qc.smk"


samples = [x.name.strip("_raw.h5ad") for x in Path(f"{PROJECT_ROOT}/raw_adata").iterdir()]

rule all:
    input: expand(f"{PROJECT_ROOT}/clean_adata/{{sample}}_denoised.h5", sample=samples)