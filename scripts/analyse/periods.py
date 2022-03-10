import xarray as xr


def periods(path_to_multiplicative_contribution_factors: str, periods: dict, path_to_output: str):
    factors = xr.open_dataset(path_to_multiplicative_contribution_factors)
    period_ds = xr.concat(
        [read_country(factors.sel(country=country), country_periods).expand_dims(country=[country])
         for country, country_periods in periods.items()],
         dim="country"
    )
    period_ds.to_netcdf(path_to_output)


def read_country(factors: xr.Dataset, periods: dict):
    return xr.concat(
        [factors.sel(year=years).mean("year").expand_dims(period=[period_name])
         for period_name, years in periods.items()],
        dim="period"
    )


if __name__ == "__main__":
    periods(
        path_to_multiplicative_contribution_factors=snakemake.input.data,
        periods=snakemake.params.periods,
        path_to_output=snakemake.output[0]
    )
