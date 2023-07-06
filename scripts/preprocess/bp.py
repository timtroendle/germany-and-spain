import pandas as pd


def bp(path_to_data: str, variables: dict[str, str], scale: bool):
    return pd.concat(
        [
            (
                pd
                .read_excel(
                    path_to_data,
                    sheet_name=variable_sheet,
                    index_col=0,
                    skiprows=2,
                    na_values=["^", "-"]
                )
                .loc[["Germany", "Spain"], 1997:2022]
                .stack()
                .pipe(scale_to_unit, sheet_name=variable_sheet, scale=scale)
                .rename_axis(index=["country", "year"])
                .rename(variable_name)
            )
            for variable_name, variable_sheet in variables.items()
        ],
        axis=1
    ).fillna(0)


def scale_to_unit(series: pd.Series, sheet_name: str, scale: bool):
    if not scale:
        return series
    if "EJ" in sheet_name:
        return series
    elif "PJ" in sheet_name:
        return series.div(1000)
    else:
        raise ValueError(f"Unknown unit of sheet {sheet_name}.")


if __name__ == "__main__":
    df = bp(
        path_to_data=snakemake.input.data,
        variables=snakemake.params.variables,
        scale=snakemake.params.scale
    )
    df.to_csv(snakemake.output[0])
