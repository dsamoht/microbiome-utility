#!/usr/bin/env python3

"""
Usage:
python kreport_mpa_barplot.py -i [KREPORT_MPA] -l [p/c/o/f/g/s] -n [INT]
python kreport_mpa_barplot.py -i [KREPORT_MPA] -l [p/c/o/f/g/s] -n [INT] -p [PARENT_TAXA_STRING]

Examples:
Plot top 15 phyla across all communities:
 -> python kreport_mpa_barplot.py -i kraken_mpa.tsv -l p -n 15

Plot top 5 orders within Cyanobacteria:
 -> python kreport_mpa_barplot.py -i kraken_mpa.tsv -l o -n 5 -p p__Cyanobacteria
"""
import argparse
import sys

from matplotlib import pyplot as plt
import pandas as pd
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as spdist

__author__ = "deschenes.thomas@gmail.com"

DISTINCT_COLORS = ['#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231',
                   '#911eb4', '#42d4f4', '#f032e6', '#bfef45', '#fabed4',
                   '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000',
                   '#aaffc3', '#808000', '#ffd8b1', '#000075', '#a9a9a9',
                   '#000000']

def format_taxa_name(lineage_str, target_level, parent_str=None):
    """
    Cleans up the Kraken MPA lineage string for the legend.
    Italicizes Genus ('g') and Species ('s') using Matplotlib's math text rendering.
    """
    if lineage_str == "Others":
        return "Others"

    nodes = lineage_str.split("|")
    target_node = next((n for n in nodes if n.startswith(f"{target_level}__")), None)

    if not target_node:
        return lineage_str

    target_name = target_node.split("__")[1].replace("_", " ")

    # Italicize Genus and Species conventions
    if target_level in ['g', 's']:
        target_name = f"$\\mathit{{{target_name}}}$"

    if parent_str:
        parent_level, parent_name = parent_str.split("__")
        parent_name = parent_name.replace("_", " ")
        if parent_level in ['g', 's']:
            parent_name = f"$\\mathit{{{parent_name}}}$"
        return f"{parent_name} | {target_name}"

    return target_name


def plot_stacked_barplot(tax_level_table, title, top_n, taxa_level, taxa_table, parent_name=None):
    """
    Shared plotting logic for clustering and generating the stacked barplot.
    """
    fig, (ax_1, ax_2) = plt.subplots(2, 1, figsize=(16, 8), gridspec_kw={'height_ratios': [1, 0.3], 'hspace': 0.6})

    # Hierarchical clustering of samples based on Bray-Curtis distance
    dist_matrix = spdist.pdist(taxa_table.T, metric='braycurtis')
    linkage_matrix = sch.linkage(dist_matrix, method='complete', optimal_ordering=True)

    dendro = sch.dendrogram(linkage_matrix,
                            orientation='bottom',
                            ax=ax_2,
                            no_labels=True,
                            color_threshold=False,
                            link_color_func=lambda x:'k')

    # Get the top_n rows by sum
    top_n_by_sum = tax_level_table.loc[tax_level_table.sum(axis=1).nlargest(top_n).index]

    # Calculate 'Others'
    if parent_name is None:
        # Relative to all
        cellular = taxa_table.loc['x__cellular_organisms'] if 'x__cellular_organisms' in taxa_table.index else 0
        viruses = taxa_table.loc['d__Viruses'] if 'd__Viruses' in taxa_table.index else 0
        others_row = (cellular + viruses) - top_n_by_sum.sum()
    else:
        # Relative to parent
        others_row = tax_level_table.loc[~tax_level_table.index.isin(top_n_by_sum.index)].sum()

    result_df = pd.concat([top_n_by_sum, pd.DataFrame(others_row).T.rename(index={0: 'Others'})])

    rel_abundance_table = result_df / result_df.sum()

    rel_abundance_table.index = [format_taxa_name(i, taxa_level, parent_name) for i in rel_abundance_table.index]

    rel_abundance_table = rel_abundance_table.T.iloc[dendro['leaves']]

    rel_abundance_table.plot(kind='bar', stacked=True, legend=False, color=DISTINCT_COLORS, ax=ax_1, edgecolor='k', width=0.8)

    n_samples = len(rel_abundance_table.index)
    ax_1.set_xlim(-0.5, n_samples - 0.5)
    ax_2.set_xlim(0, 10 * n_samples)

    ax_2.set_axis_off()
    ax_1.set_ylabel('Relative abundance')
    ax_1.legend(loc='center left', bbox_to_anchor=(1, 0.5), title='')

    fig.suptitle(title)
    fig.subplots_adjust(right=0.4)
    plt.tight_layout()
    plt.show()


def main():
    TAX_LEVELS_MAP = {'p': 6, 'c': 5, 'o': 4, 'f': 3, 'g': 2, 's': 1}
    PREFIXES = ["p__", "c__", "o__", "f__", "g__", "s__"]

    parser = argparse.ArgumentParser(description="Produce stacked barplots of microbial communities from a Kraken MPA report.")
    parser.add_argument("-i", "--input", type=str, required=True, help="Input Kraken MPA (MetaPhlAn format) file.")
    parser.add_argument("-l", "--level", type=str, required=True, choices=TAX_LEVELS_MAP.keys(), help="Target taxonomic level (p/c/o/f/g/s).")
    parser.add_argument("-n", "--top", type=int, default=10, help="Top N taxa to display in colors (max 20).")
    parser.add_argument("-p", "--parent", type=str, default=None, help="Optional exact name of parent taxa (e.g., p__Cyanobacteria).")

    args = parser.parse_args()

    if args.top > 20:
        parser.error("--top must be maximum 20 due to color palette limitations.")

    if args.parent:
        if not any(args.parent.startswith(prefix) for prefix in PREFIXES):
            parser.error("--parent must start with a valid prefix: p__, c__, o__, f__, g__, s__.")
        if TAX_LEVELS_MAP[args.parent.split("__")[0]] <= TAX_LEVELS_MAP[args.level]:
            parser.error(f"--level ({args.level}) must be a lower rank than --parent ({args.parent.split('__')[0]}).")

    try:
        main_df = pd.read_csv(args.input, sep="\t", index_col=0)
    except FileNotFoundError:
        print(f"Error: Could not find file {args.input}")
        sys.exit(1)

    if args.parent:
        tax_level_table = main_df.loc[[str(i) for i in main_df.index if args.parent in str(i) and str(i).split("|")[-1].startswith(args.level)]]
        title = f"Relative abundance of top {args.top} {args.level}-level taxa in {args.parent.split('__')[1]}"
        plot_stacked_barplot(tax_level_table, title, args.top, args.level, main_df, parent_name=args.parent)
    else:
        tax_level_table = main_df.loc[[str(i) for i in main_df.index if str(i).split("|")[-1].startswith(f"{args.level}__")]]
        title = f"Relative abundance of top {args.top} taxa at level '{args.level}'"
        plot_stacked_barplot(tax_level_table, title, args.top, args.level, main_df, parent_name=None)


if __name__ == "__main__":
    main()
