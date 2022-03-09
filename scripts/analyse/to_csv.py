import xarray as xr


def to_csv(path_to_nc, path_to_csv):
    (
        xr
        .open_dataset(path_to_nc)
        .to_dataframe()
        .to_csv(path_to_csv, header=True, index=True)
    )


if __name__ == "__main__":
    to_csv(
        path_to_nc=snakemake.input.nc,
        path_to_csv=snakemake.output[0]
    )
