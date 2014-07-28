#!/usr/bin/python3
# -*- coding: utf-8 -*-

#  Copyright 2014 Luca Beltrame <luca.beltrame@marionegri.it>
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.


import argparse
import datetime

from pathlib import Path
import joblib

from matplotlib import use

use("Qt4Agg")

from matplotlib import rc
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
import pysam
import sarge

COMMAND = "bedtools coverage -abam {0} -b {1} -counts"


pd.options.display.mpl_style = "default"
plt.ioff()

rc('font', **{'family': 'sans-serif',
              'sans-serif': ['Liberation Sans']})


def extract_samplename(bam_file):
    # Assumes only one sample per file
    samfile = pysam.Samfile(bam_file)
    return samfile.header["RG"][0]["SM"]


def calculate_coverage(sample, bed_data):

    sample_name = extract_samplename(str(sample))

    coverage = sarge.shell_format(COMMAND, sample, bed_data)

    coverage = sarge.capture_stdout(coverage)
    coverage_raw = pd.read_table(coverage.stdout, header=None,
                                 names=["chr", "start", "end", "name",
                                        "coverage"])

    # FIXME: No idea why this happens (misplaced index) but we must fix it
    if pd.Series(["+", "-"]).isin(coverage_raw["name"]).any():
        coverage_raw.reset_index(inplace=True)
        coverage_raw.columns = ["chr", "start", "end", "name", "score",
                                "strand", "coverage"]

    # clean up bed results (Split at "_", get first result)
    # FIXME: Be more flexible
    coverage_raw["name"] = coverage_raw["name"].apply(
        lambda x: x.split("_")[0].strip())

    coverage_table = coverage_raw[["coverage", "name"]].groupby("name").agg(
        ["mean", "std"])

    mean_coverage = coverage_table["coverage", "mean"]
    mean_std = coverage_table["coverage", "std"]

    coverage_table.columns.set_names(["", ""], inplace=True)
    coverage_table.index.set_names(["Gene"], inplace=True)

    coverage_table["coverage", "mean"] = mean_coverage.round(3)
    coverage_table["coverage", "std"] = mean_std.round(3)

    global_coverage_mean = coverage_raw["coverage"].mean()
    global_coverage_std = coverage_raw["coverage"].std()

    return dict(name=sample_name, cov=coverage_table,
                mean=global_coverage_mean,
                std=global_coverage_std)


def plot_coverage(coverage_data, sample_name=None):

    fig = plt.figure(figsize=(15, 10))
    axis = fig.add_subplot(1, 1, 1)
    result = coverage_data["coverage"]["mean"].plot(kind="bar", ax=axis)
    axis.set_title("Coverage data for {}".format(sample_name))

    return fig


def generate_report(coverages, tables, filename="Coverage_report.xlsx"):

    summary_df = pd.DataFrame(index=sorted(coverages.keys()),
                              columns=["Mean coverage", "Std_deviation"])
    summary_df.index.name = "Sample name"

    with pd.io.excel.ExcelWriter(filename) as excel:

        for sample in sorted(coverages):
            mean, std = coverages[sample]
            summary_df.loc[sample, "Mean coverage"] = round(mean, 4)
            summary_df.loc[sample, "Std_deviation"] = round(std, 4)

        summary_df.to_excel(excel, na_rep="NA", sheet_name="Summary")

        for sample in sorted(tables):
            dataframe = tables[sample]
            dataframe.to_excel(excel, na_rep="NA", sheet_name=sample)


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--report", action="store_true",
                        help="Write coverage report")
    parser.add_argument("source", help="Source directory containing BAM files")
    parser.add_argument("bed", help="BED file to calculate coverage for")
    parser.add_argument("destination", help="Destination file")

    options = parser.parse_args()

    source = Path(options.source)
    destination = Path(options.destination)

    if not source.exists():
        parser.error("Source directory must exist")

    files = source.glob("**/*.bam")

    coverage_calculator = joblib.delayed(calculate_coverage)
    pool = joblib.Parallel(n_jobs=-2)

    result = pool(coverage_calculator(item, options.bed)
                  for item in files)

    figures = dict()
    coverages = dict()
    tables = dict()

    for item in result:
        samplename = item["name"]
        coverage_data = item["cov"]
        figures[samplename] = plot_coverage(coverage_data, samplename)
        coverages[samplename] = (item["mean"], item["std"])
        tables[samplename] = coverage_data

    if not figures:
        print("No files found.")
        exit(1)

    if options.report:
        generate_report(coverages, tables)

    del coverages
    del tables

    with PdfPages(options.destination) as pdf:
        for figure in sorted(figures):
            pdf.savefig(figures[figure])

        d = pdf.infodict()
        d['Title'] = 'Coverage plot'
        d['Author'] = u'Parallel Genomics'
        d['Subject'] = 'Coverage plot for {} samples'.format(len(figures))
        d['CreationDate'] = datetime.datetime.today()
        d['ModDate'] = datetime.datetime.today()

if __name__ == '__main__':
    main()
