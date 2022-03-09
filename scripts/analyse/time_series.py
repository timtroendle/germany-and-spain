import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import seaborn as sns


COLOR_PALETTE = [ # Nature colors
    "#E64B35",
    "#4DBBD5",
    "#00A087",
    "#3C5488",
    "#F39B7F",
    "#8491B4",
    "#91D1C2",
    "#DC0000",
    "#7E6148",
    "#B09C85",
    "#00A087",
]
sns.set_palette(COLOR_PALETTE)


def time_series_plot(path_to_contribution_factors: str, country: str, sector: str, path_to_plot: str):
    factors = (
        xr
        .open_dataset(path_to_contribution_factors)
        .sel(country=country)
        .to_dataframe()
    )

    fig = plt.figure(figsize=(8, 3))
    ax = fig.add_subplot(111)

    factors.plot(ax=ax)

    ax.set_ylabel("Change since crisis")
    ax.get_xaxis().set_major_locator(MultipleLocator(10))
    ax.get_xaxis().set_minor_locator(MultipleLocator(1))
    ax.set_title(f"{country} -- {sector}")
    sns.despine(fig, right=False)
    fig.tight_layout()
    fig.savefig(path_to_plot)


if __name__ == "__main__":
    time_series_plot(
        path_to_contribution_factors=snakemake.input.data,
        country=snakemake.wildcards.country,
        sector=snakemake.wildcards.sector,
        path_to_plot=snakemake.output[0]
    )
