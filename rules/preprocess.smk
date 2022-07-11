rule extract_data_from_excel:
    message: "Extract all input data from Excel file."
    input:
        data = config["data-sources"]["excel"]["path"]
    params:
        dataset_config = config["data-sources"]["excel"]["datasets"]
    output: "build/raw-data.nc"
    conda: "../envs/preprocess.yaml"
    script: "../scripts/preprocess/extract.py"


rule factors:
    message: "Derive factors from raw data."
    input:
        data = rules.extract_data_from_excel.output[0]
    output:
        total = "build/factors-total.nc",
        industry = "build/factors-industry.nc",
        transport = "build/factors-transport.nc",
        power = "build/factors-power.nc"
    conda: "../envs/preprocess.yaml"
    script: "../scripts/preprocess/factors.py"


rule download_primary_energy:
    message: "Download primary energy data."
    params: url = config["data-sources"]["primary-energy"]["url"]
    output: protected("data/automatic/raw-primary-energy.csv")
    shell: "curl -sLo {output} '{params.url}'"


rule primary_energy:
    message: "Preprocess primary energy data."
    input:
        data = rules.download_primary_energy.output[0]
    params:
        countries = COUNTRIES,
        years = range(config["parameters"]["first-year"], config["parameters"]["final-year"] + 1),
        columns = config["data-sources"]["primary-energy"]["columns"]
    output: "build/primary-energy-twh.csv"
    conda: "../envs/preprocess.yaml"
    script: "../scripts/preprocess/primary_energy.py"
