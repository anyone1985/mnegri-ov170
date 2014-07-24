#!/usr/bin/env python
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

from __future__ import print_function

import argparse
from itertools import groupby
import os
import sys

from gemini import GeminiQuery
from pathlib import Path
import matplotlib
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import pandas as pd
import sarge
import cruzdb

# This uses the public UCSC instance. In case MySQL is blocked, you can
# point it using a valid SQLalchemy URL to an existing database server.

ucsc = cruzdb.Genome(db="hg19")

pd.options.display.mpl_style = "default"
rc('font', **{'family': 'sans-serif',
              'sans-serif': ['Liberation Sans']})

renames = {"ts": "transition", "tv": "transversion", "ins": "insertion",
           "del": "deletion"}


def generate_phenotypes(database):

    query = GeminiQuery(database)
    query_string = "SELECT name, phenotype FROM samples"
    phenotypes = {1: list(), 2: list()}

    query.run(query_string)

    for row in query:
        phenotypes[int(row["phenotype"])].append(row["name"])

    return phenotypes


def get_nearby_gene(chrom, start, end):

    nearest = ucsc.knearest("refFlat", chrom, start, end)
    assert len(nearest) == 1
    nearest = nearest[0].geneName

    return nearest


def check_multiple_alts(chrom, start, end,  alt, samples, reference_file=None):

    # FIXME: PyVCF can't be used as it loads the wrong coordinate with
    # fetch()

    if reference_file is None:
        return alt

    alt = alt.split(",")

    if len(alt) == 1:
        return alt[0]

    region = "{0}:{1}-{2}".format(chrom, start, end)
    sample_display = ",".join(samples)

    bcftools = ("/usr/bin/bcftools view {0} -r {1} -a -s {2}"
                " --exclude-uncalled -H")
    bcftools = sarge.shell_format(bcftools, reference_file, region,
                                  sample_display)
    command = sarge.capture_both(bcftools)

    mutation_table = pd.read_table(command.stdout, header=None,
                                   names=["CHROM", "POS", "ID", "REF", "ALT"],
                                   usecols=["CHROM", "POS", "ID", "REF",
                                            "ALT"])

    if mutation_table.empty:
        # Try alternative approach: sometimes bcftools asserts when
        # launched from sarge
        import subprocess
        with open(os.devnull) as null:
            cmd = subprocess.Popen(bcftools, shell=True,
                                   stdout=subprocess.PIPE,
                                   stderr=null)
        mutation_table = pd.read_table(cmd.stdout, header=None,
                                       names=["CHROM", "POS", "ID",
                                              "REF", "ALT"],
                                       usecols=["CHROM", "POS", "ID",
                                                "REF", "ALT"])

    seen_alt = mutation_table["ALT"].item()

    if len(seen_alt.split(",")) > 1:
        message = ("Detected more than one allele for sample pair {}: {},"
                   " mutation position {}-{}, ref {}")
        print(message.format(", ".join(samples), seen_alt, chrom, start,
                             mutation_table["REF"].item()),
              file=sys.stderr)
        # This is a false positive! Same locus but different ALTs
        return

    return seen_alt


def extract_shared_mutations(database, reference_file=None):

    phenotypes = generate_phenotypes(database)
    query = GeminiQuery(database)
    query_string = ("SELECT chrom, start, end, gene, ref, alt, type, sub_type,"
                    "impact, codon_change, aa_change, vcf_id, cosmic_ids"
                    " FROM variants WHERE in_1kg=0")

    query.run(query_string, show_variant_samples=True)
    rows = list()

    for row in query:

        variants = row.variant_samples

        if any(item in phenotypes[1] for item in variants) and any(
               item in phenotypes[2] for item in variants):

            valid_groups = list()

            chrom, start, end, alt = (row["chrom"], row["start"],
                                      row["end"], row["alt"])

            # In the case of intergenic regions, get the name of the
            # closest gene

            if row["gene"] is None:
                gene = get_nearby_gene(chrom, start, end)
                # print "None subsituted with", gene
            else:
                gene = row["gene"]

            for gid, group in groupby(variants, lambda x: x.split("_")[1]):

                # Rename according to Pandora guidelines
                # Starts with 0: 1 + number
                # 3 digits: 20 + number
                # 4 digits: 2 + number

                if len(gid) < 5:
                    newgid = "2" + gid if len(gid) == 4 else "20" + gid
                else:
                    newgid = gid

                group = list(group)

                if len(list(group)) == 2:

                    # Check if we have different ALT bases for the samples
                    # in the same pair. If this occurs, it is a false positive
                    # and should be discarded. To do so, we need a VCF file to
                    # query by base, otherwise we take the value as-is.

                    if reference_file is not None:

                        alt = check_multiple_alts(chrom, start, end, alt,
                                                  group, reference_file)

                        if alt is None:
                            # Biallelic site for pair - discard
                            continue

                    valid_groups.append(newgid)

            cosmic_data = "Yes" if row["cosmic_ids"] else "No"

            data = [chrom, start, end, gene, row["ref"], alt, row["type"],
                    row["sub_type"], row["impact"], row["codon_change"],
                    row["aa_change"], row["vcf_id"], cosmic_data]

            if not valid_groups:
                rows.append(data + [np.nan])
            else:
                for gid in valid_groups:
                    rows.append(data + [gid])

    colnames = ["chrom", "start", "end", "gene", "ref", "alt", "type",
                "sub_type", "impact", "codon_change", "aa_change",
                "dbsnp_id", "in_cosmic", "variants_with_pairs"]

    df = pd.DataFrame.from_records(rows, columns=colnames)
    df.set_index(["chrom", "start", "end", "gene", "ref", "alt"], inplace=True)
    # Get rid of loci without pairs
    df = df.dropna(subset=["variants_with_pairs"])

    return df


def colorbar_text(ax, text):

    plt.setp(ax.get_yticklabels(), visible=False)
    plt.setp(ax.get_xticklabels(), visible=False)

    ax.axis('off')
    ax.grid(False)

    ax.text(1.2, 0.25, text, horizontalalignment="right")


def _set_axis_parameters(ax, dataframe, names=True):

    # put the major ticks at the middle of each cell
    ax.set_yticks(np.arange(dataframe.shape[0]) + 0.5,
                  minor=False)
    ax.set_xticks(np.arange(dataframe.shape[1]) + 0.5,
                  minor=False)
    ax.set_xlim((0, dataframe.shape[1]))
    ax.set_ylim((0, dataframe.shape[0]))

    ax.invert_yaxis()
    ax.xaxis.tick_top()

    plt.xticks(rotation=90)

    if names:
        ax.set_xticklabels(dataframe.columns, minor=False, rotation=90)
    else:
        ax.set_xticklabels([])

    ax.set_yticklabels(dataframe.index, minor=False)

    # Turn off all the ticks
    ax = plt.gca()

    for t in ax.xaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    for t in ax.yaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False

    ax.grid(False)


def clean_axis(ax):

    """Remove ticks, tick labels, and frame from axis"""
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)


def colorbar_text(ax, text):

    ax.axis('off')
    ax.grid(False)

    ax.text(1.2, 0.25, text, horizontalalignment="right")


def generate_heatmap(dataframe_with_pairs, samplelist=None,
                     pathway_colortable=None,
                     histotype_colortable=None):

    import matplotlib.gridspec as grid

    plt.ioff()

    dataframe = dataframe_with_pairs.reset_index()
    grouper = dataframe.groupby("gene")
    # Get counts per gene and per sample
    counter = grouper["variants_with_pairs"].value_counts()

    # index.levels[0] is the gene, index.levels[1] are pairs
    # The new structure is a table with genes as rows and
    # sample pairs as columns
    new_dataframe = pd.DataFrame(index=counter.index.levels[0],
                                 columns=counter.index.levels[1])
    new_dataframe.fillna(0, inplace=True)

    # Put the counts in their proper place: we iterate over a MultiIndex of
    # gene, pair couples, along with the count (the value), and we just
    # use the first two as a coordinate where to put the third

    for rowid, row in counter.iteritems():
        new_dataframe.loc[rowid[0], rowid[1]] = row

    if samplelist is not None:

        # Put missing samples in
        missing = samplelist[~samplelist.isin(new_dataframe.columns)]
        for item in missing:
            new_dataframe[item] = 0

        new_dataframe.sort_index(axis=1, inplace=True)

    # Add pathway elements

    fig = plt.figure(figsize=(15, 10))

    if pathway_colortable is not None and histotype_colortable is not None:

        present_genes = pathway_colortable[
            pathway_colortable.index.isin(new_dataframe.index)]
        present_samples = histotype_colortable.loc[new_dataframe.index]

        # Same order as the color table
        new_dataframe = new_dataframe[histotype_colortable.index.astype("str")]

        gs = grid.GridSpec(3, 2, width_ratios=[0.5, 15],
                           height_ratios=[0.2, 0.2, 15])

        pathway_ax = fig.add_subplot(gs[2, 0])  # Pathways
        heatmap_ax = fig.add_subplot(gs[2, 1], sharey=pathway_ax)  # Heatmap

        histotype_ax = fig.add_subplot(gs[0, 1])  # Histotype bar
        grade_ax = fig.add_subplot(gs[1, 1])  # Grade bar

        # Those stay in the top row, and represent the labels next
        # to the color bars

        bartext_ax = fig.add_subplot(gs[0, 0])  # Histotype
        gradetext_ax = fig.add_subplot(gs[1, 0])  # Grade

        clean_axis(histotype_ax)
        clean_axis(grade_ax)
        clean_axis(bartext_ax)
        clean_axis(gradetext_ax)

        # These will be our labels on top of histotype and grade
        create_colorbar(histotype_colortable, histotype_ax, False,
                        labels=True)
        create_colorbar(histotype_colortable, grade_ax, False, "Grade_color",
                        labels=False)

        # Labels for histotype and grade

        colorbar_text(bartext_ax, "Histotype")
        colorbar_text(gradetext_ax, "Grade")

        new_dataframe = new_dataframe.loc[present_genes.index]
        plt.setp(heatmap_ax.get_yticklabels(), visible=False)
        create_colorbar(present_genes, pathway_ax)

    else:
        fig.add_subplot(1, 1, 1)
        heatmap_ax = fig.axes[0]

    cmap, norm = create_colormap()

    heatmap = heatmap_ax.pcolor(new_dataframe, cmap=cmap,
                                edgecolors="black", alpha=1, norm=norm)

    labels = True if pathway_colortable is None else False

    _set_axis_parameters(heatmap_ax, new_dataframe, names=labels)

    cax = fig.add_axes([-0.05, 1.025, 0.15, 0.025])

    cbar = fig.colorbar(heatmap, cax=cax, orientation="horizontal",
                        ticks=range(5))  # HACK: Hardcoded!
    cbar.solids.set_edgecolor("face")

    if pathway_colortable is not None:
        plt.tight_layout()

    return fig, cax


def _reduce_value(value):

    if value.name != "samples":
        # Get the unique name, should be identical
        return value.unique().item()

    return ", ".join(value)


def _correct_names(name):

    name = name.split("_")[1]

    if len(name) < 5:
        name = "2" + name if len(name) == 4 else "20" + name

    return name


def create_colorbar(dataframe, ax, row=True, color_col="Color", labels=True):

    gene_names = dataframe.index
    colors = dataframe[color_col]
    length = list(range(len(dataframe.index)))
    matrix, cmap = color_list_to_matrix_and_cmap(colors, length, row=row)
    mesh = ax.pcolormesh(matrix, cmap=cmap)

    if not row:
        ax.set_xlim(0, matrix.shape[1])
        ax.set_ylim((0, matrix.shape[0]))
        ax.set_xticks(np.arange(matrix.shape[1]) + 0.5,
                      minor=False)
    else:
        ax.set_xlim(0, matrix.shape[1])
        ax.set_ylim((0, matrix.shape[0]))
        ax.set_yticks(np.arange(matrix.shape[0]) + 0.5,
                      minor=False)

    if not row:
        ax.invert_yaxis()
        ax.xaxis.tick_top()

        if labels:
            ax.set_xticklabels(gene_names, minor=False, rotation=90)
        else:
            ax.set_xticklabels([])

        ax.set_yticks([])
    else:

        # We invert gene_names because in this case the x axis is flipped
        # (heatmap)
        if labels:
            ax.set_yticklabels(gene_names[::-1], minor=False)
        else:
            ax.set_yticklabels([])
        ax.set_xticks([])

    return mesh


def color_list_to_matrix_and_cmap(colors, ind, row=True):

    """Turns a list of colors into a numpy matrix and matplotlib colormap
    For 'heatmap()'
    This only works for 1-column color lists..

    These arguments can now be plotted using matplotlib.pcolormesh(matrix,
    cmap) and the provided colors will be plotted.

    Parameters
    ----------
    colors : list of matplotlib colors
        Colors to label the rows or columns of a dataframe.
    ind : list of ints
        Ordering of the rows or columns, to reorder the original colors
        by the clustered dendrogram order
    row : bool
        Is this to label the rows or columns? Default True.

    Returns
    -------
    matrix : numpy.array
        A numpy array of integer values, where each corresponds to a color
        from the originally provided list of colors
    cmap : matplotlib.colors.ListedColormap

    """
    # TODO: Support multiple color labels on an element in the heatmap
    import matplotlib as mpl

    colors_original = colors
    colors = set(colors)
    col_to_value = dict((col, i) for i, col in enumerate(colors))
    matrix = np.array([col_to_value[col] for col in colors_original])[ind]

    # Is this row-side or column side?
    if row:
        # shape of matrix: nrows x 1
        new_shape = (len(colors_original), 1)
    else:
        # shape of matrix: 1 x ncols
        new_shape = (1, len(colors_original))
    matrix = matrix.reshape(new_shape)

    cmap = mpl.colors.ListedColormap(colors)
    return matrix, cmap


def create_colormap():

    # HACK: max is hardcoded!
    norm = colors.Normalize(vmin=0, vmax=4)
    # define the colormap
    cmap = plt.cm.autumn
    # extract all colors from the .jet map
    cmaplist = [cmap(i) for i in range(cmap.N)]
    # force the first color entry to be grey
    cmaplist[0] = (.5, .5, .5, 1.0)
    # create the new map
    cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)

    return cmap, norm


def save_figure(fig, destination, extra_artist=None):

    from pathlib import Path
    name = Path(destination)

    if extra_artist is not None:
        extra_args = {"extra_artists": (extra_artist, ),
                      "bbox_inches": "tight"}
    else:
        extra_args = {}

    for extension in [".pdf", ".svgz", ".png"]:
        if extension != "png":
            fig.savefig(str(name.with_suffix(extension)), **extra_args)
        else:
            fig.savefig(str(name.with_suffix(extension)), dpi=300,
                        **extra_args)


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--reference-file", help="Reference VCF file")
    parser.add_argument("--heatmap", help="Plot a heatmap",
                        action="store_true")
    parser.add_argument("-c", "--color-table",
                        help="Add a pathway color table")
    parser.add_argument("-f", "--phenotype-table",
                        help="Add a phenotype/grade color table")
    parser.add_argument("--keep-synonymous", action="store_true",
                        help="Keep synonymous mutations in the table")
    parser.add_argument("database", help="Source GEMINI database")
    parser.add_argument("destination", help="Destination table to save to")

    options = parser.parse_args()
    destination = Path(options.destination)

    print("Extracting shared mutations...", file=sys.stderr)
    df = extract_shared_mutations(options.database, options.reference_file)

    if not options.keep_synonymous:
        df = df[df.impact != "synonymous_coding"]
    else:
        print("Including synonymous mutations.", file=sys.stderr)

    df.to_csv(str(destination), sep="\t", na_rep="NA")
    print("Saved text file: {}".format(destination), file=sys.stderr)

    # User-visible file: rename the column to make more sense
    xls_df = df.rename(columns={"variants_with_pairs": "samples"})

    grouper = xls_df.groupby(level=[0, 1, 2, 3, 4, 5])

    # Put all common samples for one mutation on one line
    xls_df = grouper.aggregate(_reduce_value)
    # Rename for clearer understanding
    xls_df.sub_type = xls_df.sub_type.apply(
        lambda x: renames.get(x, "unknown"))

    xlsx_destination = destination.with_suffix(".xlsx")
    xls_df.to_excel(str(xlsx_destination), sheet_name="Result",
                    na_rep="NA", merge_cells=False)
    print("Saved xlsx file: {}".format(xlsx_destination), file=sys.stderr)

    if options.heatmap:

        samplelist = generate_phenotypes(options.database)
        samplelist = pd.Series(samplelist[1])
        samplelist = samplelist.apply(_correct_names).drop_duplicates()

        if options.color_table is not None:
            colortable = pd.read_table(options.color_table, index_col=0)
        else:
            colortable = None

        if options.phenotype_table is not None:
            phenotype_table = pd.read_table(options.phenotype_table,
                                            index_col=0)
        else:
            phenotype_table = None

        fig, cax = generate_heatmap(df, samplelist, colortable,
                                    phenotype_table)
        save_figure(fig, destination, extra_artist=cax)

if __name__ == '__main__':
    main()
