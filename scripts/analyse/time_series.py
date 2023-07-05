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
GREY = "#7F7F7F"
CRISES = ["global-financial", "covid"]


def time_series_plot_all_countries(path_to_contribution_factors: str, path_to_plot: str,
                                   periods: dict[str, dict[str, list[int]]]):
    ds = (
        xr
        .open_dataset(path_to_contribution_factors)
    )

    fig = plt.figure(figsize=(8, 4))
    axes = fig.subplots(len(ds.country), 1, sharex=True, sharey=True)
    for i, (ax, country) in enumerate(zip(axes, ds.country)):
        country = country.item()
        df = (
            ds
            .sel(country=country, drop=True)
            .to_dataframe()
            .rename(columns=lambda name: name.replace("-", " ").capitalize().replace("Gdp", "GDP"))
            .rename_axis(index="Year")
        )
        time_series_plot_single_countries(
            factors=df,
            country=country,
            ax=ax,
            legend=True if i == 0 else False,
            periods=periods[country]
        )

    sns.despine(fig)
    fig.tight_layout()
    fig.savefig(path_to_plot)


def time_series_plot_single_countries(factors: pd.DataFrame, country: str, ax: plt.Axes, legend: bool,
                                      periods: dict[str, list[int]]):
    ax.set_prop_cycle(color=COLOR_PALETTE[:6], linestyle=['-', '--', '-.', ':', '--', '-.'])
    factors["Emissions"].plot(ax=ax, linewidth=2.5, color=COLOR_PALETTE[0])
    factors.drop(columns=["Emissions"]).plot(ax=ax)
    ax.set_ylabel("Change since 1998")
    ax.get_xaxis().set_major_locator(MultipleLocator(2))
    ax.get_xaxis().set_minor_locator(MultipleLocator(1))
    ax.set_title(f"{country}")

    for i, crisis in enumerate(CRISES):
        ax.axvspan(
            xmin=periods[crisis][0],
            xmax=periods[crisis][-1] + 1,
            ymin=0,
            ymax=1,
            linewidth=0.0,
            alpha=0.2,
            color=GREY,
            label="Crises" if i == 0 else None
        )

    if legend:
        ax.legend(bbox_to_anchor=(1, 1), loc="upper left", frameon=False)
    else:
        ax.legend().remove()


if __name__ == "__main__":
    time_series_plot_all_countries(
        path_to_contribution_factors=snakemake.input.data,
        periods=snakemake.params.periods,
        path_to_plot=snakemake.output[0]
    )
