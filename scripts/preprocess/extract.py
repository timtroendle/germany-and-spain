import pandas as pd
import xarray as xr


def extract(path_to_data: str, dataset_config: list[list], path_to_output: str):
    variables = [
        read_var(path_to_data, sheet_name, name, first_row)
        for sheet_name, name, first_row in dataset_config
    ]
    ds = xr.Dataset({var.name: var for var in variables})
    ds.to_netcdf(path_to_output)



def read_var(path_to_file: str, sheet_name: str, name: str, first_row: int):
    return (
        pd
        .read_excel(path_to_file, sheet_name=sheet_name, skiprows=first_row - 1, nrows=2, index_col=0)
        .rename_axis(columns="year", index="country")
        .transpose()
        .unstack()
        .to_xarray()
        .rename(f"{sheet_name}-{name}")
    )


if __name__ == "__main__":
    extract(
        path_to_data=snakemake.input.data,
        dataset_config=snakemake.params.dataset_config,
        path_to_output=snakemake.output[0]
    )
