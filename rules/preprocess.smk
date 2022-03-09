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
