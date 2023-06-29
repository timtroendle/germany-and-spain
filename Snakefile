from snakemake.utils import min_version

PANDOC = "pandoc --filter pantable --filter pandoc-fignos --filter pandoc-tablenos --citeproc"

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
    message: "Run entire analysis and compile report."
    input:
        expand(
            "build/figures/{country}-{sector}-time-series.png",
            country=COUNTRIES,
            sector=SECTORS
        ),
        expand(
            "build/figures/{country}-sectoral-emissions.png",
            country=COUNTRIES
        ),
        expand(
            "build/figures/{country}-primary-energy.png",
            country=COUNTRIES
        ),
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


def pandoc_options(wildcards):
    suffix = wildcards["suffix"]
    if suffix == "html":
        return "--self-contained --to html5"
    elif suffix == "pdf":
        return "--pdf-engine weasyprint"
    elif suffix == "docx":
        return []
    else:
        raise ValueError(f"Cannot create report with suffix {suffix}.")


rule report:
    message: "Compile report.{wildcards.suffix}."
    input:
        "report/literature.yaml",
        "report/report.md",
        "report/pandoc-metadata.yaml",
        "report/apa.csl",
        "report/reset.css",
        "report/report.css",
    params: options = pandoc_options
    output: "build/report.{suffix}"
    wildcard_constraints:
        suffix = "((html)|(pdf)|(docx))"
    conda: "envs/report.yaml"
    shadow: "minimal"
    shell:
        """
        cd report
        ln -s ../build .
        {PANDOC} report.md  --metadata-file=pandoc-metadata.yaml {params.options} \
        -o ../build/report.{wildcards.suffix}
        """


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
