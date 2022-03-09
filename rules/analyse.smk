rule multiplicative_contribution:
    message: "Calculate multiplicative contribution of all factors in sector {wildcards.sector}."
    input:
        script = "scripts/analyse/contribution.py",
        data = "build/factors-{sector}.nc"
    output: "build/multiplicative-contribution-factors-{sector}.nc"
    conda: "../envs/default.yaml"
    script: "../scripts/analyse/contribution.py"


rule relative_cumulative_contributions:
    message: "Calculate relative cumulative contributions of all factors in sector {wildcards.sector}."
    input:
        script = "scripts/analyse/cumulative.py",
        data = "build/multiplicative-contribution-factors-{sector}.nc"
    output: "build/relative-cumulative-contribution-factors-{sector}.nc"
    conda: "../envs/default.yaml"
    script: "../scripts/analyse/cumulative.py"


rule contribution_time_series_plot:
    message: "Plot contribution time series of all factors in sector {wildcards.sector} "
             + "in country {wildcards.country}."
    input:
        script = "scripts/analyse/time_series.py",
        data = "build/relative-cumulative-contribution-factors-{sector}.nc"
    output: "build/{country}-{sector}-time-series.png"
    conda: "../envs/default.yaml"
    script: "../scripts/analyse/time_series.py"
