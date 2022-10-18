rule multiplicative_contribution:
    message: "Calculate multiplicative contribution of all factors in sector {wildcards.sector}."
    input:
        data = "build/factors-{sector}.nc"
    output: "build/multiplicative-contribution-factors-{sector}.nc"
    conda: "../envs/default.yaml"
    script: "../scripts/analyse/contribution.py"


rule relative_cumulative_contributions:
    message: "Calculate relative cumulative contributions of all factors in sector {wildcards.sector}."
    input:
        data = "build/multiplicative-contribution-factors-{sector}.nc"
    output: "build/relative-cumulative-contribution-factors-{sector}.nc"
    conda: "../envs/default.yaml"
    script: "../scripts/analyse/cumulative.py"


rule periods:
    message: "Calculate cumulative contribution during crises in sector {wildcards.sector}."
    input:
        data = "build/multiplicative-contribution-factors-{sector}.nc"
    params:
        periods = config["parameters"]["periods"]
    output: "build/periods-{sector}.nc"
    conda: "../envs/default.yaml"
    script: "../scripts/analyse/periods.py"


rule contribution_time_series_plot:
    message: "Plot contribution time series of all factors in sector {wildcards.sector} "
             + "in country {wildcards.country}."
    input:
        data = "build/relative-cumulative-contribution-factors-{sector}.nc"
    output: "build/{country}-{sector}-time-series.png"
    conda: "../envs/default.yaml"
    script: "../scripts/analyse/time_series.py"


rule primary_energy_plot:
    message: "Plot sources of primary energy over time in country {wildcards.country}."
    input:
        data = rules.primary_energy.output[0]
    output: "build/{country}-primary-energy.png"
    conda: "../envs/default.yaml"
    script: "../scripts/analyse/primary_energy.py"


rule emissions_plot:
    message: "Plot sectoral emission time series."
    input:
        industry = "build/factors-industry.nc",
        transport = "build/factors-transport.nc",
        power = "build/factors-power.nc"
    output: "build/{country}-sectoral-emissions.png"
    wildcard_constraints:
        country = "|".join(COUNTRIES)
    conda: "../envs/default.yaml"
    script: "../scripts/analyse/sectoral_emissions.py"


rule netcdf_to_csv:
    message: "Transform {wildcards.filename}.nc to csv."
    input:
        nc = "build/{filename}.nc"
    output: "build/{filename}.csv"
    conda: "../envs/default.yaml"
    script: "../scripts/analyse/to_csv.py"
