import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import seaborn as sns


COLOURS = {
    "Oil": "#5d5d5d",
    "Gas": "#b9b9b9",
    "Coal": "#181818",
    "Nuclear": "#cc0000",
    "Biofuel": "#8fce00",
    "Hydro": "#2986cc",
    "Solar": "#ffd966",
    "Wind": "#674ea7",
    "Other renewable": "#e062db",
}


def primary_energy_all_countries(path_to_data: str, path_to_plot: str):
    ds = (
        pd
        .read_csv(path_to_data, index_col=[0, 1])
        .to_xarray()
    )

    fig = plt.figure(figsize=(8, 4))
    axes = fig.subplots(len(ds.country), 1, sharex=True, sharey=True)
    for i, (ax, country) in enumerate(zip(axes, ds.country)):
        df = (
            ds
            .sel(country=country, drop=True)
            .to_dataframe()
            .rename(columns=str.capitalize)
        )
        primary_energy_single_country(df, country.item(), ax, legend=True if i == 0 else False)

    sns.despine(fig)
    fig.tight_layout()
    fig.savefig(path_to_plot)


def primary_energy_single_country(df: pd.DataFrame, country: str, ax: plt.Axes, legend: bool):
    ax.stackplot(
        df.index.values,
        df.to_dict(orient="list").values(),
        colors=[COLOURS[tech] for tech in df.columns.values],
        labels=df.columns.values,
        alpha=0.7
    )
    if legend:
        ax.legend(bbox_to_anchor=(1, 1), loc="upper left", frameon=False)
    ax.set_title(country)
    ax.set_ylabel("Primary energy\nconsumption (EJ)")
    ax.get_xaxis().set_major_locator(MultipleLocator(2))
    ax.get_xaxis().set_minor_locator(MultipleLocator(1))

    ax.get_yaxis().set_major_locator(MultipleLocator(5))
    ax.get_yaxis().set_minor_locator(MultipleLocator(1))
    if ax.get_ylim()[1] < 10:
        ax.set_ylim(top=10)


if __name__ == "__main__":
    primary_energy_all_countries(
        path_to_data=snakemake.input.data,
        path_to_plot=snakemake.output[0]
    )
