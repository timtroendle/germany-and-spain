import eurostat
import pandas as pd
import pycountry


def download_eurostat_data(code: str, name: str, filter_pars: dict[str, list[str]],
                           cat_mapping: dict[str, dict[str: str]], col_name_mapping: dict[str, str],
                           additional_index_cols: list[str]) -> pd.DataFrame:
    return (
        eurostat
        .get_data_df(code, filter_pars=filter_pars)
        .assign(**{
                column: lambda df: df[column].map(mapping)
                for column, mapping in cat_mapping.items()
        })
        .assign(
            country=lambda df: df["geo\TIME_PERIOD"].map(lambda c: pycountry.countries.lookup(c).name)
        )
        .rename(columns=col_name_mapping)
        .set_index(["country"] + additional_index_cols)
        .filter(axis="columns", regex="\d\d\d\d")
        .rename_axis(columns="year")
        .stack()
        .rename(name)
    )


if __name__ == "__main__":
    df = download_eurostat_data(
        code=snakemake.params["eurostat-dataset-code"],
        name=snakemake.params["name"],
        filter_pars=snakemake.params["filter-pars"],
        additional_index_cols=snakemake.params["additional-index-cols"],
        cat_mapping=snakemake.params["cat-mapping"],
        col_name_mapping=snakemake.params["col-name-mapping"]
    )
    df.to_csv(snakemake.output[0], index=True, header=True)
