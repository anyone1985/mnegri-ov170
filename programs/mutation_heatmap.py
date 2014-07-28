#!/usr/bin/python3
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
from pathlib import Path

import fastcluster
import numpy as np
import numpy.ma as ma
import pandas as pd

from matplotlib import rc
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as grid
from matplotlib.colors import Normalize
import scipy.cluster.hierarchy as sch


pd.options.display.mpl_style = "default"
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Liberation Sans']})

# These are needed due to amplicons spanning over the allotted area.

substitutions = {
    "ZEB1-AS1": "ZEB1",      # Different strand
    "RHNO1": "FOXM1",        # Different strand
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


def save_figure(fig: plt.Figure, destination: str,
                extra_artist: mpl.artist.Artist=None):

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

    This function was taken from https://github.com/olgabot/seaborn and
    is under a BSD license.

    """
    # TODO: Support multiple color labels on an element in the heatmap

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


def _set_colorbar_axis_limits(ax: mpl.axis.Axis,
                              matrix: pd.DataFrame,
                              row: bool,
                              index: pd.Series,
                              labels: bool=True):

    if not row:
        ax.set_xlim(0, matrix.shape[1])
        ax.set_ylim((0, matrix.shape[0]))
        ax.set_xticks(np.arange(matrix.shape[1]) + 0.5,
                      minor=False)
        ax.invert_yaxis()
        ax.xaxis.tick_top()

        # Column mode: the gene_names list is in the right order
        if labels:
            ax.set_xticklabels(index, minor=False, rotation=90)
        else:
            ax.set_xticklabels([])

        ax.set_yticks([])

    else:
        ax.set_xlim(0, matrix.shape[1])
        ax.set_ylim((0, matrix.shape[0]))
        ax.set_yticks(np.arange(matrix.shape[0]) + 0.5,
                      minor=False)

        # Row mode: the heat map flips x axis, so we need to reverse the
        # colorbar
        if labels:
            ax.set_yticklabels(index[::-1], minor=False)
        else:
            ax.set_yticklabels([])
        ax.set_xticks([])


def create_colorbar(dataframe: pd.DataFrame,
                    ax: mpl.axis.Axis, row: bool=True,
                    color_col: str="Color",
                    labels: bool=True) -> mpl.artist.Artist:

    index = dataframe.index

    continuous = False
    norm = None
    matrix_mask = None
    shading = None

    if isinstance(color_col, str):
        colors = dataframe[color_col]
        length = list(range(len(dataframe.index)))
        matrix, cmap = color_list_to_matrix_and_cmap(colors, length, row=row)

    elif isinstance(color_col, list) or isinstance(color_col, tuple):

        matrix, cmap = color_col
        matrix = matrix.loc[dataframe.index]

        norm = Normalize(vmin=matrix.min(), vmax=matrix.max())

        if row:
            matrix = matrix.reshape(len(matrix), 1)
        else:
            matrix = matrix.reshape(1, len(matrix))

        # Missing values need masked arrays!
        if pd.isnull(matrix).any():
            matrix_mask = ma.array(matrix, mask=np.isnan(matrix))

        continuous = True

    _set_colorbar_axis_limits(ax, matrix, row, index, labels)

    if continuous and matrix_mask is not None:

        # left, right, top, bottom
        extent = [ax.get_xlim()[0], ax.get_xlim()[1],
                  ax.get_ylim()[1], ax.get_ylim()[0]]

        # if matrix_mask is not None:
        cmap = cm.get_cmap(cmap)
        cmap.set_bad("grey")
        matrix = matrix_mask

    else:
        mesh = ax.pcolormesh(matrix, cmap=cmap, norm=norm)
        mesh.set_edgecolor("face")

    return mesh


def colorbar_text(ax: mpl.axis.Axis, text: str):

    plt.setp(ax.get_yticklabels(), visible=False)
    plt.setp(ax.get_xticklabels(), visible=False)

    ax.axis('off')
    ax.grid(False)

    ax.text(1.4, 0.25, text, horizontalalignment="right")


def create_colormap(max=8):

    norm = Normalize(vmin=0, vmax=max)
    # define the colormap
    cmap = plt.cm.autumn
    # extract all colors from the .jet map
    cmaplist = [cmap(i) for i in range(cmap.N)]
    # force the first color entry to be grey
    cmaplist[0] = (.5, .5, .5, 1.0)
    # create the new map
    cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)

    return cmap, norm


def set_axis_text(ax: mpl.axis.Axis, dataframe: pd.DataFrame):

    # put the major ticks at the middle of each cell
    ax.set_yticks(np.arange(dataframe.shape[0]) + 0.5,
                  minor=False)
    ax.set_xticks(np.arange(dataframe.shape[1]) + 0.5,
                  minor=False)
    ax.set_xlim((0, dataframe.shape[1]))
    ax.set_ylim((0, dataframe.shape[0]))


def set_axis_parameters(ax: mpl.axis.Axis, dataframe: pd.DataFrame,
                        names: bool=True):

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


def clean_axis(ax: mpl.axis.Axis):

    """Remove ticks, tick labels, and frame from axis"""
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)


def generate_heatmap(dataframe: pd.DataFrame,
                     color_dataframe: pd.DataFrame,
                     histo_dataframe: pd.DataFrame,
                     color_map: dict=None,
                     cluster: bool=True,
                     figsize: tuple=(10, 15)) -> plt.Figure:

    fig = plt.figure(figsize=figsize)

    if cluster:
        linkage = fastcluster.linkage(dataframe.T,
                                      "complete",
                                      metric="correlation",
                                      preserve_input=True)
        dendrogram_row_ratio = 2
    else:
        linkage = None
        # Make row smaller without dendrogram
        dendrogram_row_ratio = 0.1

    max_rows = 4 if color_map is not None else 3

    # The bar plot is actually on the third row
    # FIXME: Handle long labels (ratios)

    if color_map is not None:
        bar_ratio = 0.25 * len(color_map)
        height_ratios = [dendrogram_row_ratio, 0.25, bar_ratio, 15]
    else:
        height_ratios = [dendrogram_row_ratio, 0.25, 15]

    gs = grid.GridSpec(max_rows, 2, height_ratios=height_ratios,
                       width_ratios=[0.2, 15])

    dendro_ax = fig.add_subplot(gs[0, 1], axisbg="white")  # Dendrogram

    plt.setp(dendro_ax.get_yticklabels(), visible=False)

    pathway_ax = fig.add_subplot(gs[-1, 0])  # Pathway
    heatmap_ax = fig.add_subplot(gs[-1, 1], sharey=pathway_ax)  # Heatmap

    # Con gridspec e' necessario fare questo in modo che gli assi Y
    # non siano visibili

    plt.setp(heatmap_ax.get_yticklabels(), visible=False)

    if linkage is not None:
        leaf_ax = fig.add_subplot(gs[1, 1], sharex=dendro_ax)
        dendro = sch.dendrogram(linkage, ax=dendro_ax, no_labels=False,
                                labels=dataframe.columns,
                                leaf_rotation=90,)
        # Reorder dataframe according to the labels in the leaves
        dataframe = dataframe[dendro["ivl"]]  # Leaf node labels

        # Put labels in the right order!
        histo_dataframe = histo_dataframe.loc[dendro["ivl"]]

        # TRICK: Given that printing labels screws layout because they add an
        # x axis, we generate a specific axis only with the text, iterating on
        # the locations of the labels of the dendrogram. After the new text is
        # in place, we remove the labels from the dendrogram.

        for leafname, leafcoord in zip(dendro["ivl"],
                                       dendro_ax.xaxis.get_ticklocs()):
            leaf_ax.text(leafcoord, 0.99, leafname, rotation=90,
                         horizontalalignment="center")

    else:
        set_axis_parameters(heatmap_ax, dataframe, False)
        leaf_ax = fig.add_subplot(gs[1, 1], sharex=heatmap_ax)
        dataframe = dataframe.loc[:, histo_dataframe.index]

        set_axis_parameters(leaf_ax, dataframe, False)

        for index, leafcoord in enumerate(leaf_ax.xaxis.get_ticklocs()):
            leaf_ax.text(leafcoord, 0.99, dataframe.columns[index],
                         rotation=90, horizontalalignment="center")

    clean_axis(leaf_ax)
    leaf_ax.grid(False)
    leaf_ax.axis('off')
    clean_axis(dendro_ax)

    if color_map is not None:

        subgrids = len(color_map)
        gs_inside = grid.GridSpecFromSubplotSpec(
            subgrids, 1,
            subplot_spec=gs[2, 1],
            height_ratios=[1 for item in color_map])

        bars = list()
        for index, group in enumerate(sorted(color_map)):

            column = color_map[group]
            bar_ax = fig.add_subplot(gs_inside[index])
            clean_axis(bar_ax)
            create_colorbar(histo_dataframe, bar_ax, False, column,
                            labels=False)
            bar_ax.text(-0.25, 0.5, group, horizontalalignment="right",
                        verticalalignment="center")
            bars.append(bar_ax)

    create_colorbar(color_dataframe, pathway_ax)

    cmap, norm = create_colormap()

    dataframe = dataframe.loc[color_dataframe.index]

    heatmap1 = heatmap_ax.pcolor(dataframe, cmap=cmap,
                                 edgecolors="black", alpha=1,
                                 norm=norm)

    set_axis_parameters(heatmap_ax, dataframe, False)

    cax = fig.add_axes([-0.05, 1.025, 0.15, 0.025])

    cbar = fig.colorbar(heatmap1, cax=cax, orientation="horizontal",
                        ticks=range(9))

    cbar.solids.set_edgecolor("face")

    gs.tight_layout(fig)

    return fig, cax


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--color-table",
                        help="Add a pathway color table")
    parser.add_argument("-f", "--phenotype-table",
                        help="Add a phenotype/grade color table")
    parser.add_argument("-t", "--type", help="Sample type (fts or sds)",
                        choices=("fts", "sds"), default="fts")
    parser.add_argument("--cluster", action="store_true",
                        help="Perform hierarchical clustering")
    parser.add_argument("source", help="Source mutation table")
    parser.add_argument("destination", help="File to save the plot to")

    options = parser.parse_args()

    frame = pd.read_table(options.source, index_col=0)

    try:
        frame.drop("PLA2G6", inplace=True)  # Unwanted
    except ValueError:
        pass

    frame.rename(substitutions, inplace=True)
    frame = frame.groupby(level=0).sum()  # Remove duplicates and sum

    if options.color_table is not None:
        colortable = pd.read_table(options.color_table, index_col=0)
    else:
        colortable = None

    if options.phenotype_table is not None:
        phenotype_table = pd.read_table(options.phenotype_table,
                                        index_col=0)
    else:
        phenotype_table = None

    # Map colors to the right column

    color_map = {"Histotype": "Color",
                 "Grade": "Grade_color",
                 "PFS": "PFS_state",
                 "OS": "OS_state"}

    if options.type == "fts":
        color_map["Pt-sensitivity"] = "Resistance_color_fts"
    else:
        color_map["Pt-sensitivity"] = "Resistance_color"

    fig, cax = generate_heatmap(frame, options.color_table,
                                options.phenotype_table,
                                color_map=color_map, cluster=options.cluster,
                                figsize=(15, 15))

    save_figure(fig, options.destination, extra_artist=cax)

if __name__ == '__main__':
    main()
