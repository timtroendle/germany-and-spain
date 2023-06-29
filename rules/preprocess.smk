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


rule extract_primary_energy_from_excel:
    message: "Extract primary energy data from Excel file."
    input:
        data = config["data-sources"]["primary-energy"]["path"]
    params:
        dataset_config = config["data-sources"]["primary-energy"]["datasets"]
    output: "build/data/raw-primary-energy.nc"
    conda: "../envs/preprocess.yaml"
    script: "../scripts/preprocess/extract.py"


rule primary_energy:
    message: "Preprocess primary energy data."
    input:
        data = rules.extract_primary_energy_from_excel.output[0]
    output: "build/data/primary-energy-ej.csv"
    conda: "../envs/preprocess.yaml"
    script: "../scripts/preprocess/primary_energy.py"
