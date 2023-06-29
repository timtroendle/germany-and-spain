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

FIRST_YEAR = 1998


def time_series_plot_all_countries(path_to_industry: str, path_to_power: str, path_to_transport: str,
                                   path_to_plot: str):
    industry, power, transport = [
        (
            xr
            .open_dataset(path)["emissions"]
        )
        for path in [path_to_industry, path_to_power, path_to_transport]
    ]
    fig = plt.figure(figsize=(8, 4))
    axes = fig.subplots(len(industry.country), 1, sharex=True, sharey=True)
    for i, (ax, country) in enumerate(zip(axes, industry.country)):
        time_series = [
            (
                ds
                .sel(country=country, drop=True)
                .to_dataframe()
                .rename_axis(index="Year")
                .loc[FIRST_YEAR:]
            )
            for ds in [industry, power, transport]
        ]
        time_series_plot_single_country(*time_series, country.item(), ax, legend=True if i == 0 else False)

    sns.despine(fig)
    fig.tight_layout()
    fig.savefig(path_to_plot)


def time_series_plot_single_country(industry: pd.DataFrame, power: pd.DataFrame, transport: pd.DataFrame,
                                    country: str, ax: plt.Axes, legend: bool):
    ax.set_prop_cycle(color=COLOR_PALETTE[:6], linestyle=['-', '--', '-.', ':', '--', '-.'])
    ax.plot(industry.div(industry.iloc[0]), label="Industry")
    ax.plot(power.div(power.iloc[0]), label="Power")
    ax.plot(transport.div(transport.iloc[0]), label="Transport")
    ax.legend()

    ax.set_ylabel(f"Change since {FIRST_YEAR}")
    ax.get_xaxis().set_major_locator(MultipleLocator(10))
    ax.get_xaxis().set_minor_locator(MultipleLocator(1))
    ax.set_title(country)
    if legend:
        ax.legend(bbox_to_anchor=(1, 1), loc="upper left", frameon=False)
    else:
        ax.legend().remove()


if __name__ == "__main__":
    time_series_plot_all_countries(
        path_to_industry=snakemake.input.industry,
        path_to_power=snakemake.input.power,
        path_to_transport=snakemake.input.transport,
        path_to_plot=snakemake.output[0]
    )
