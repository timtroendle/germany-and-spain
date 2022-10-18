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


def time_series_plot(path_to_industry: str, path_to_power: str, path_to_transport: str,
                     country: str, path_to_plot: str):
    industry, power, transport = [
        (
            xr
            .open_dataset(path)["emissions"]
            .sel(country=country)
            .to_dataframe()["emissions"]
            .rename_axis(index="Year")
            .loc[FIRST_YEAR:]
        )
        for path in [path_to_industry, path_to_power, path_to_transport]
    ]

    fig = plt.figure(figsize=(8, 3))
    ax = fig.add_subplot(111)

    ax.set_prop_cycle(color=COLOR_PALETTE[:6], linestyle=['-', '--', '-.', ':', '--', '-.'])
    ax.plot(industry.div(industry.iloc[0]), label="Industry")
    ax.plot(power.div(power.iloc[0]), label="Power")
    ax.plot(transport.div(transport.iloc[0]), label="Transport")
    ax.legend()

    ax.set_ylabel(f"Change since {FIRST_YEAR}")
    ax.get_xaxis().set_major_locator(MultipleLocator(10))
    ax.get_xaxis().set_minor_locator(MultipleLocator(1))
    ax.set_title(f"Sectoral emissions in {country}")
    sns.despine(fig)
    fig.tight_layout()
    fig.savefig(path_to_plot)


if __name__ == "__main__":
    time_series_plot(
        path_to_industry=snakemake.input.industry,
        path_to_power=snakemake.input.power,
        path_to_transport=snakemake.input.transport,
        country=snakemake.wildcards.country,
        path_to_plot=snakemake.output[0]
    )
