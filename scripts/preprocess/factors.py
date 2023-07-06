import pandas as pd
import xarray as xr


def factors(emissions: xr.DataArray, population: xr.DataArray, gdp: xr.DataArray, energy: xr.DataArray) -> xr.Dataset:
    return (
        xr
        .Dataset({
            "emissions": emissions,
            "population": population,
            "affluence": gdp / population,
            "energy-intensity": energy / gdp,
            "carbon-intensity": emissions / energy
        })
    )


def read_file(path_to_file: str) -> xr.DataArray:
    return (
        pd
        .read_csv(path_to_file, index_col=[0, 1])
        .to_xarray()
        .to_array()
        .isel(variable=0, drop=True)
    )


if __name__ == "__main__":
    ds = factors(
        emissions=read_file(snakemake.input.emissions),
        gdp=read_file(snakemake.input.gdp),
        energy=read_file(snakemake.input.energy),
        population=read_file(snakemake.input.population),
    )
    ds.to_netcdf(snakemake.output[0])
