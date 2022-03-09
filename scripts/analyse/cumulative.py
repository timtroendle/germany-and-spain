import xarray as xr


def cumulative_factors(path_to_factors: str, path_to_output: str):
    factors = (
        xr
        .open_dataset(path_to_factors)
        .dropna("year")
    )
    cumulative = factors.cumprod("year")
    cumulative["year"] = factors["year"]
    rel_factors = cumulative / cumulative.isel(year=0)
    rel_factors.to_netcdf(path_to_output)


if __name__ == "__main__":
    cumulative_factors(
        path_to_factors=snakemake.input.data,
        path_to_output=snakemake.output[0]
    )
