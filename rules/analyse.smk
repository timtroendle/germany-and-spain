rule factors:
    message: "Derive factors from raw data."
    input:
        emissions = rules.emissions.output[0],
        energy = rules.primary_energy.output[0],
        population = rules.population.output[0],
        gdp = rules.gdp.output[0]
    output: "build/results/factors.nc"
    conda: "../envs/preprocess.yaml"
    script: "../scripts/preprocess/factors.py"


rule multiplicative_contribution:
    message: "Calculate multiplicative contribution of all factors."
    input:
        data = "build/results/factors.nc"
    output: "build/results/multiplicative-contribution-factors.nc"
    conda: "../envs/default.yaml"
    script: "../scripts/analyse/contribution.py"


rule relative_cumulative_contributions:
    message: "Calculate relative cumulative contributions of all factors."
    input:
        data = "build/results/multiplicative-contribution-factors.nc"
    params:
        first_year = config["parameters"]["first-year"]
    output: "build/results/relative-cumulative-contribution-factors.nc"
    conda: "../envs/default.yaml"
    script: "../scripts/analyse/cumulative.py"


rule periods:
    message: "Calculate cumulative contribution during crises."
    input:
        data = "build/results/multiplicative-contribution-factors.nc"
    params:
        periods = config["parameters"]["periods"]
    output: "build/results/periods.nc"
    conda: "../envs/default.yaml"
    script: "../scripts/analyse/periods.py"



rule contribution_time_series_plot:
    message: "Plot contribution time series of all factors."
    input:
        data = "build/results/relative-cumulative-contribution-factors.nc"
    params: periods = config["parameters"]["periods"]
    output: "build/figures/contribution-time-series.png"
    conda: "../envs/default.yaml"
    script: "../scripts/analyse/time_series.py"


rule primary_energy_by_fuel_plot:
    message: "Plot primary energy by fuel over time."
    input:
        data = rules.primary_energy_by_fuel.output[0]
    output: "build/figures/primary-energy-by-fuel.png"
    conda: "../envs/default.yaml"
    script: "../scripts/analyse/primary_energy.py"


rule sectoral_emissions_plot:
    message: "Plot sectoral emission time series."
    input:
        emissions = rules.sectoral_emissions.output[0]
    params:
        first_year = config["parameters"]["first-year"],
        periods = config["parameters"]["periods"]
    output: "build/figures/sectoral-emissions.png"
    conda: "../envs/default.yaml"
    script: "../scripts/analyse/sectoral_emissions.py"


rule netcdf_to_csv:
    message: "Transform {wildcards.filename}.nc to csv."
    input:
        nc = "build/{pathname}/{filename}.nc"
    output: "build/{pathname}/{filename}.csv"
    conda: "../envs/default.yaml"
    script: "../scripts/analyse/to_csv.py"
