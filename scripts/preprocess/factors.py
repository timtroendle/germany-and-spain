import xarray as xr


def factors(path_to_data: str, paths_to_output: list[str]):
    ds = xr.open_dataset(path_to_data)
    (
        xr
        .Dataset({
            "emissions": ds["Total-emissions"],
            "population": ds["Total-population"],
            "affluence": ds["Total-gdp"] / ds["Total-population"],
            "energy-intensity": ds["Total-energy"] / ds["Total-gdp"],
            "carbon-intensity": ds["Total-emissions"] / ds["Total-energy"]
        })
        .to_netcdf(paths_to_output.total)
    )
    (
        xr
        .Dataset({
            "emissions": ds["Industry-emissions"],
            "gdp": ds["Total-gdp"],
            "gdp-share": ds["Industry-gdp"] / ds["Total-gdp"],
            "energy-intensity": ds["Industry-energy"] / ds["Industry-gdp"],
            "carbon-intensity": ds["Industry-emissions"] / ds["Industry-energy"]
        })
        .to_netcdf(paths_to_output.industry)
    )
    (
        xr
        .Dataset({
            "emissions": ds["Transport-emissions"],
            "population": ds["Total-population"],
            "affluence": ds["Total-gdp"] / ds["Total-population"],
            "energy-intensity": ds["Transport-energy"] / ds["Total-gdp"],
            "carbon-intensity": ds["Transport-emissions"] / ds["Transport-energy"]
        })
        .to_netcdf(paths_to_output.transport)
    )
    (
        xr
        .Dataset({
            "emissions": ds["Power-emissions"],
            "population": ds["Total-population"],
            "affluence": ds["Total-gdp"] / ds["Total-population"],
            "efficiency": ds["Power-energy-in"] / ds["Power-energy-out"],
            "energy-intensity": ds["Power-energy-out"] / ds["Total-gdp"],
            "carbon-intensity": ds["Power-emissions"] / ds["Power-energy-in"]
        })
        .to_netcdf(paths_to_output.power)
    )


if __name__ == "__main__":
    factors(
        path_to_data=snakemake.input.data,
        paths_to_output=snakemake.output
    )
