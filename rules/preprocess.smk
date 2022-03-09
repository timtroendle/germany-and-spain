rule extract_data_from_excel:
    message: "Extract all input data from Excel file."
    input:
        script = "scripts/preprocess/extract.py",
        data = config["data-sources"]["excel"]["path"]
    params:
        dataset_config = config["data-sources"]["excel"]["datasets"]
    output: "build/raw-data.nc"
    conda: "../envs/preprocess.yaml"
    script: "../scripts/preprocess/extract.py"


rule factors:
    message: "Derive factors from raw data."
    input:
        script = "scripts/preprocess/factors.py",
        data = rules.extract_data_from_excel.output[0]
    output:
        total = "build/factors-total.nc",
        industry = "build/factors-industry.nc",
        transport = "build/factors-transport.nc",
        power = "build/factors-power.nc"
    conda: "../envs/preprocess.yaml"
    script: "../scripts/preprocess/factors.py"
