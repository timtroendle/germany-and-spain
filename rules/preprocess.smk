rule download_bp_data:
    message: "Download BP data file."
    params: url = config["data-sources"]["bp"]["url"]
    output:
        protected("data/automatic/raw-bp.xslx")
    shell:
        "curl -sLo {output} '{params.url}'"


rule primary_energy:
    message: "Preprocess primary energy data."
    input:
        data = rules.download_bp_data.output[0]
    params:
        variables = config["data-sources"]["primary-energy"]["variables"],
        scale = config["data-sources"]["primary-energy"]["scale"]
    output: "build/data/primary-energy-ej.csv"
    conda: "../envs/preprocess.yaml"
    script: "../scripts/preprocess/bp.py"


rule primary_energy_by_fuel:
    message: "Preprocess primary energy data by fuel."
    input:
        data = rules.download_bp_data.output[0]
    params:
        variables = config["data-sources"]["primary-energy-by-fuel"]["variables"],
        scale = config["data-sources"]["primary-energy-by-fuel"]["scale"]
    output: "build/data/primary-energy-by-fuel-ej.csv"
    conda: "../envs/preprocess.yaml"
    script: "../scripts/preprocess/bp.py"


rule emissions:
    message: "Preprocess emission data."
    input:
        data = rules.download_bp_data.output[0]
    params:
        variables = config["data-sources"]["emissions"]["variables"],
        scale = config["data-sources"]["emissions"]["scale"]
    output: "build/data/emissions-mt.csv"
    conda: "../envs/preprocess.yaml"
    script: "../scripts/preprocess/bp.py"


rule sectoral_emissions:
    message: "Download and preprocess sectoral emissions data."
    params:
        **config["data-sources"]["sectoral-emissions"]
    output: protected("build/data/raw-sectoral-emissions-mt.csv")
    conda: "../envs/preprocess.yaml"
    script: "../scripts/preprocess/eurostat.py"


rule gdp:
    message: "Download and preprocess GDP."
    params:
        **config["data-sources"]["gdp"]
    output: protected("build/data/raw-gdp-meur.csv")
    conda: "../envs/preprocess.yaml"
    script: "../scripts/preprocess/eurostat.py"


rule population:
    message: "Download and preprocess population."
    params:
        **config["data-sources"]["population"]
    output: protected("build/data/raw-population.csv")
    conda: "../envs/preprocess.yaml"
    script: "../scripts/preprocess/eurostat.py"
