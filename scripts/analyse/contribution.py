import xarray as xr


def contribution(path_to_factors: str, path_to_output: str):
    factors = xr.open_dataset(path_to_factors)
    multiplicative_contribution = factors / factors.shift({"year": 1})
    multiplicative_contribution.to_netcdf(path_to_output)


if __name__ == "__main__":
    contribution(
        path_to_factors=snakemake.input.data,
        path_to_output=snakemake.output[0]
    )
