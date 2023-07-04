import pandas as pd


def primary_energy(path_to_data: str, carriers: dict[str, str]):
    return pd.concat(
        [
            (
                pd
                .read_excel(
                    path_to_data,
                    sheet_name=carrier_sheet,
                    index_col=0,
                    skiprows=2,
                    na_values=["^", "-"]
                )
                .loc[["Germany", "Spain"], 1997:2022]
                .stack()
                .pipe(scale_to_unit, carrier_sheet)
                .rename_axis(index=["country", "year"])
                .rename(carrier_name)
            )
            for carrier_name, carrier_sheet in carriers.items()
        ],
        axis=1
    ).fillna(0)


def scale_to_unit(series: pd.Series, sheet_name: str):
    if "EJ" in sheet_name:
        return series
    elif "PJ" in sheet_name:
        return series.div(1000)
    else:
        raise ValueError(f"Unknown unit of sheet {sheet_name}.")


if __name__ == "__main__":
    df = primary_energy(
        path_to_data=snakemake.input.data,
        carriers=snakemake.params.carriers
    )
    df.to_csv(snakemake.output[0])
