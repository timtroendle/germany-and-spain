from snakemake.utils import min_version

COUNTRIES = ["Germany", "Spain"]
SECTORS = ["total", "industry", "transport", "power"]

configfile: "config/default.yaml"
include: "rules/preprocess.smk"
include: "rules/analyse.smk"
min_version("7.8")


onsuccess:
    if "email" in config.keys():
        shell("echo "" | mail -s 'germany-and-spain succeeded' {config[email]}")
onerror:
    if "email" in config.keys():
        shell("echo "" | mail -s 'germany-and-spain failed' {config[email]}")

rule all:
    message: "Run all analysis steps."
    input:
        expand(
            "build/figures/{sector}-time-series.png",
            country=COUNTRIES,
            sector=SECTORS
        ),
        "build/figures/sectoral-emissions.png",
        "build/figures/primary-energy.png",
        expand(
            "build/results/relative-cumulative-contribution-factors-{sector}.csv",
            sector=SECTORS
        ),
        expand(
            "build/results/multiplicative-contribution-factors-{sector}.csv",
            sector=SECTORS
        ),
        expand(
            "build/results/periods-{sector}.csv",
            sector=SECTORS
        )


rule dag:
     message: "Plot dependency graph of the workflow."
     output:
         dot = "build/dag.dot",
         pdf = "build/dag.pdf"
     conda: "envs/dag.yaml"
     shell:
         """
         snakemake --rulegraph > {output.dot}
         dot -Tpdf -o {output.pdf} {output.dot}
         """


rule push:
    message: "Package, zip, and move entire build."
    params: push_directory = config["push-directory"]
    shell:
        """
        zip -r {params.push_directory}/germany-and-spain-$(date -Idate).zip build
        """


rule clean: # removes all generated results
    message: "Remove all build results but keep downloaded data."
    run:
         import shutil

         shutil.rmtree("build")
         print("Data downloaded to data/ has not been cleaned.")


rule test:
    conda: "envs/test.yaml"
    output: "build/test-report.html"
    shell:
        "py.test --html={output} --self-contained-html"
