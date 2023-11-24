import pandas as pd
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
MARKER_SIZE = 3


def time_series_plot_all_countries(emissions: pd.DataFrame, path_to_plot: str, first_year: int,
                                   periods: dict[str, dict[str, list[int]]]):
    emissions = emissions.to_xarray()
    fig = plt.figure(figsize=(8, 4))
    axes = fig.subplots(len(emissions.country), 1, sharex=True, sharey=True)
    for i, (ax, country) in enumerate(zip(axes, emissions.country)):
        country = country.item()
        normed_time_series = (
            emissions
            .sel(country=country, year=slice(1998, None), drop=True)
        )
        normed_time_series = normed_time_series / normed_time_series.isel(year=0)

        time_series_plot_single_country(
            normed_time_series.to_dataframe()["emissions_mt"],
            country=country,
            ax=ax,
            legend=True if i == 0 else False,
            periods=periods[country],
            first_year=first_year
        )

    sns.despine(fig)
    fig.tight_layout()
    fig.savefig(path_to_plot)


def time_series_plot_single_country(df: pd.DataFrame, country: str, ax: plt.Axes, legend: bool,
                                    first_year: int, periods: dict[str, list[int]]):
    ax.set_prop_cycle(color=COLOR_PALETTE[:6], linestyle=['-', '--', '-.', ':', '--', '-.'])
    (
        df
        .unstack(0)
        .rename(columns=lambda name: name.capitalize())
        .rename_axis(index="Year")
        .plot(ax=ax, marker="o", markersize=MARKER_SIZE)
    )

    ax.set_ylabel(f"Change since {first_year}")
    ax.get_xaxis().set_major_locator(MultipleLocator(2))
    ax.get_xaxis().set_minor_locator(MultipleLocator(1))
    ax.set_title(country)

    for i, crisis in enumerate(CRISES):
        ax.axvspan(
            xmin=periods[crisis][0] - 0.5,
            xmax=periods[crisis][-1] + 0.5,
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
        emissions=pd.read_csv(snakemake.input.emissions, index_col=[0, 1, 2]),
        path_to_plot=snakemake.output[0],
        periods=snakemake.params.periods,
        first_year=snakemake.params.first_year
    )
