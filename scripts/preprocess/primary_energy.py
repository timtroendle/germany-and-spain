import pandas as pd


def primary_energy(path_to_data: str, countries: list[str], years: list[int],
                   columns: list[str], path_to_output: str):
    (
        pd
        .read_csv(path_to_data, index_col=[1, 2])[columns]
        .to_xarray()
        .sel(country=countries, year=years)
        .to_dataframe()
        .to_csv(path_to_output, index=True, header=True)
    )


if __name__ == "__main__":
    primary_energy(
        path_to_data=snakemake.input.data,
        countries=snakemake.params.countries,
        years=snakemake.params.years,
        columns=snakemake.params.columns,
        path_to_output=snakemake.output[0]
    )
