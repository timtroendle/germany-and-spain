rule extract_data_from_excel:
    message: "Extract all input data from Excel file."
    input:
        data = config["data-sources"]["factors"]["path"]
    params:
        dataset_config = config["data-sources"]["factors"]["datasets"]
    output: "build/data/raw.nc"
    conda: "../envs/preprocess.yaml"
    script: "../scripts/preprocess/extract.py"


rule factors:
    message: "Derive factors from raw data."
    input:
        data = rules.extract_data_from_excel.output[0]
    output:
        total = "build/data/factors-total.nc",
        industry = "build/data/factors-industry.nc",
        transport = "build/data/factors-transport.nc",
        power = "build/data/factors-power.nc"
    conda: "../envs/preprocess.yaml"
    script: "../scripts/preprocess/factors.py"


rule download_primary_energy:
    message: "Download primary energy file."
    params: url = config["data-sources"]["primary-energy"]["url"]
    output:
        protected("data/automatic/raw-primary-energy.xslx")
    shell:
        "curl -sLo {output} '{params.url}'"


rule primary_energy:
    message: "Preprocess primary energy data."
    input:
        data = rules.download_primary_energy.output[0]
    params:
        carriers = config["data-sources"]["primary-energy"]["carriers"]
    output: "build/data/primary-energy-ej.csv"
    conda: "../envs/preprocess.yaml"
    script: "../scripts/preprocess/primary_energy.py"


rule emissions:
    message: "Download and preprocess emissions data."
    params:
        **config["data-sources"]["emissions"]
    output: "data/automatic/raw-emissions-mt.csv"
    conda: "../envs/preprocess.yaml"
    script: "../scripts/preprocess/eurostat.py"
