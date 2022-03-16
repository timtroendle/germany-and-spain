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


def primary_energy(path_to_data: str, country: str, path_to_plot: str):
    df = (
        pd
        .read_csv(path_to_data, index_col=[0, 1])
        .to_xarray()
        .sel(country=country)
        .to_dataframe()
        .drop(columns=["country"])
        .rename(columns=lambda name: name.replace("_consumption", "").replace("_", " ").capitalize())
    )

    fig = plt.figure(figsize=(8, 3))
    ax = fig.add_subplot(111)

    ax.stackplot(
        df.index.values,
        df.div(1000).to_dict(orient="list").values(),
        colors=[COLOURS[tech] for tech in df.columns.values],
        labels=df.columns.values,
        alpha=0.7
    )
    ax.legend(bbox_to_anchor=(1.0, 1.0), loc="upper left")
    ax.set_title(country)
    ax.set_ylabel("Primary energy consumption (PWh)")
    ax.get_xaxis().set_major_locator(MultipleLocator(10))
    ax.get_xaxis().set_minor_locator(MultipleLocator(1))

    ax.get_yaxis().set_major_locator(MultipleLocator(1))
    ax.get_yaxis().set_minor_locator(MultipleLocator(0.25))
    if ax.get_ylim()[1] < 2:
        ax.set_ylim(top=2)

    sns.despine(fig)
    fig.tight_layout()
    fig.savefig(path_to_plot)


if __name__ == "__main__":
    primary_energy(
        path_to_data=snakemake.input.data,
        country=snakemake.wildcards.country,
        path_to_plot=snakemake.output[0]
    )
