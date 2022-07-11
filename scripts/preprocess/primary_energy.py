import xarray as xr


def primary_energy(path_to_data: str, path_to_output: str):
    (
        xr
        .open_dataset(path_to_data)
        .to_dataframe()
        .rename(columns=lambda col: col.split("-")[-1])
        .to_csv(path_to_output, index=True, header=True)
    )


if __name__ == "__main__":
    primary_energy(
        path_to_data=snakemake.input.data,
        path_to_output=snakemake.output[0]
    )
