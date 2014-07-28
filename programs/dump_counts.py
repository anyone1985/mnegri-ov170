#!/usr/bin/env python
# -*- coding: utf-8 -*-

#  Copyright 2014 Luca Beltrame <luca.beltrame@marionegri.it>
#
#  This file is part of utils.
#
#  utils is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  utils is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with utils.  If not, see <http://www.gnu.org/licenses/>.

import argparse
from collections import defaultdict, Counter

import itertools

from gemini import GeminiQuery
import cruzdb
import pandas as pd

# This uses the public UCSC instance. In case MySQL is blocked, you can
# point it using a valid SQLalchemy URL to an existing database server.

ucsc = cruzdb.Genome(db="hg19")

# These work around proper gene association bugs (due to amplicon extending
# past boundaries, or annotation picking things on the wrong strand)

substitutions = {
    "ZEB1-AS1": "ZEB1",      # Different strand
    "RHNO1": "FOXM1",        # Different strand
    "PLA2G6": np.nan,     # off-target?
    "C3orf72": "FOXL2",      # Different strand
    # "MC1R": np.nan,       # past the promoter
    "ACAA1": "MYD88",        # Different strand
    "VIM-AS1": "VIM",        # Different strand
    "LOC100507424": "FOXM1",  # Wrong annotation?
    "MTOR-AS1": "MTOR",
    "EGFR-AS1": "EGFR",
    "WRAP53": "TP53",
    "EPM2AIP1": "MLH1",
    "C5orf22": "DROSHA",
    "C9orf53": "CDKN2A",
    "LYRM5": "KRAS",
    "N4BP2L1": "BRCA2",
    "RMDN3": "RAD51",
    "NBR2": "BRCA1",
    "CNTD2": "AKT2",
    "HSCB": "CHEK2",
    "NPAT": "ATM",
    "MC1R": "TUBB3"
}


def get_nearby_gene(chrom, start, end):

    nearest = ucsc.knearest("refFlat", chrom, start, end)
    nearest = pd.Series([item.geneName for item in nearest])
    nearest = nearest.apply(lambda x: substitutions[x]
                            if x in substitutions else x)
    try:
        nearest = nearest.drop_duplicates().item()
    except Exception:
        print(nearest.drop_duplicates())
        raise

    # assert len(nearest) == 1

    return nearest


def summarize_by_gene_and_sample(db, coding_only=True):

    "This is copied from GEMINI's own burden tool"

    query = ("select chrom, start, end, gene, impact, info from variants where"
             " impact != 'synonymous_coding' and in_1kg=0 ")

    if coding_only:
        query += " and codon_change != 'None'"

    gq = GeminiQuery(db)
    gq.run(query, show_variant_samples=True)

    burden = defaultdict(Counter)

    counter = itertools.count(1)
    current_count = 0

    for row in gq:

        gene_name = row['gene']

        if not gene_name:
            gene_name = get_nearby_gene(row["chrom"], row["start"],
                                        row["end"])

        new_counts = Counter(row["HET_samples"])
        # Counter can't do scalar multiplication
        new_counts = new_counts + Counter(row["HOM_ALT_samples"])
        new_counts = new_counts + Counter(row["HOM_ALT_samples"])

        del new_counts['']
        burden[gene_name] += new_counts

    dfs = list()
    for gene_name, counts in burden.items():
        df = pd.DataFrame(counts, columns=counts.keys(),
                          index=[gene_name])
        dfs.append(df)

    df = pd.concat(dfs)
    df = df.fillna(0)

    return df


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--coding-only", action="store_true")
    parser.add_argument("source", help="Source GEMINI database")
    parser.add_argument("destination", help="destination file to save to")

    options = parser.parse_args()
    df = summarize_by_gene_and_sample(options.source, options.coding_only)

    with open(options.destination, "w") as handle:
        df.to_csv(handle, sep="\t", na_rep="NA")


if __name__ == '__main__':
        main()