import xarray as xr


def cumulative_factors(path_to_factors: str, reference_year: int, path_to_output: str):
    factors = (
        xr
        .open_dataset(path_to_factors)
        .dropna("year")
    )
    cumulative = factors.cumprod("year")
    cumulative["year"] = factors["year"]
    cumulative = cumulative.sel(year=slice(reference_year, None))
    rel_factors = cumulative / cumulative.isel(year=0)
    rel_factors.to_netcdf(path_to_output)


if __name__ == "__main__":
    cumulative_factors(
        path_to_factors=snakemake.input.data,
        reference_year=snakemake.params.first_year,
        path_to_output=snakemake.output[0]
    )
