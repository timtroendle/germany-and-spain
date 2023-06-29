rule multiplicative_contribution:
    message: "Calculate multiplicative contribution of all factors in sector {wildcards.sector}."
    input:
        data = "build/data/factors-{sector}.nc"
    output: "build/results/multiplicative-contribution-factors-{sector}.nc"
    conda: "../envs/default.yaml"
    script: "../scripts/analyse/contribution.py"


rule relative_cumulative_contributions:
    message: "Calculate relative cumulative contributions of all factors in sector {wildcards.sector}."
    input:
        data = "build/results/multiplicative-contribution-factors-{sector}.nc"
    output: "build/results/relative-cumulative-contribution-factors-{sector}.nc"
    conda: "../envs/default.yaml"
    script: "../scripts/analyse/cumulative.py"


rule periods:
    message: "Calculate cumulative contribution during crises in sector {wildcards.sector}."
    input:
        data = "build/results/multiplicative-contribution-factors-{sector}.nc"
    params:
        periods = config["parameters"]["periods"]
    output: "build/results/periods-{sector}.nc"
    conda: "../envs/default.yaml"
    script: "../scripts/analyse/periods.py"


rule contribution_time_series_plot:
    message: "Plot contribution time series of all factors in sector {wildcards.sector}."
    input:
        data = "build/results/relative-cumulative-contribution-factors-{sector}.nc"
    params: periods = config["parameters"]["periods"]
    output: "build/figures/{sector}-time-series.png"
    conda: "../envs/default.yaml"
    script: "../scripts/analyse/time_series.py"


rule primary_energy_plot:
    message: "Plot sources of primary energy over time."
    input:
        data = rules.primary_energy.output[0]
    output: "build/figures/primary-energy.png"
    conda: "../envs/default.yaml"
    script: "../scripts/analyse/primary_energy.py"


rule emissions_plot:
    message: "Plot sectoral emission time series."
    input:
        industry = "build/data/factors-industry.nc",
        transport = "build/data/factors-transport.nc",
        power = "build/data/factors-power.nc"
    params: periods = config["parameters"]["periods"]
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
