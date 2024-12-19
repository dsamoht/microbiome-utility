#!/usr/bin/env python3

"""
usage:
python kreport_mpa_barplot.py [KREPORT_MPA] [COLLAPSED_SEQTAB]
"""
import argparse

from matplotlib import pyplot as plt
import pandas as pd
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as spdist


__author__ = "deschenes.thomas@gmail.com"

# source : https://sashamaps.net/docs/resources/20-colors/
DISTINCT_COLORS = ['#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231',
                   '#911eb4', '#42d4f4', '#f032e6', '#bfef45', '#fabed4',
                   '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000',
                   '#aaffc3', '#808000', '#ffd8b1', '#000075', '#a9a9a9',
                   '#000000']

def barplot_relative_to_all(taxa_level, top_n, taxa_table):
    """
    Stacked barplot of the relative abundance of `top_n` `taxa_level`
    `taxa_table`
    
    ex: Usage for a stacked barplot of top 15 phylum:
     -> barplot_relative_to_all('p', 15, kraken_mpa_table)
    """
    
    fig, (ax_1, ax_2) = plt.subplots(2, 1, figsize=(16, 8), gridspec_kw={'height_ratios': [1, 0.3], 'hspace': 0.6})
    
    tax_level_table = taxa_table.loc[[str(i) for i in taxa_table.index if str(i).split("|")[-1].startswith(f"{taxa_level}__")]]
    
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
    others_row = (taxa_table.loc['x__cellular_organisms'] + taxa_table.loc['d__Viruses'])  - top_n_by_sum.sum()
    result_df = pd.concat([top_n_by_sum, pd.DataFrame(others_row).T.rename(index={0: 'Others'})])
    
    tax_level_table = result_df / result_df.sum()
    tax_level_table.index = ["|".join([j for j in i.split("|") if "x__" not in j]) for i in tax_level_table.index]
    tax_level_table = tax_level_table.T.iloc[dendro['leaves']]
    
    tax_level_table.plot(kind='bar', stacked=True, legend=False, color=DISTINCT_COLORS, ax=ax_1, edgecolor='k', width=0.8)    
    ax_2.set_axis_off()
    ax_2.set_line_color = 'black'
    ax_1.set_ylabel('Relative abundance')
    ax_1.legend(
            loc='center left', 
            bbox_to_anchor=(1, 0.5), 
            title=''
        )
    fig.suptitle(f"Relative abundance of top {top_n} {taxa_level}")
    fig.subplots_adjust(right=0.4)
    plt.tight_layout()
    plt.show()

def barplot_relative_to_parent(taxon_name, taxa_level, top_n, taxa_table):
    
    """
    Stacked barplot of `top_n` `division` inside `taxon_name`
    
    ex: Usage for a stacked barplot of top 5 Cyanobacterial orders:
     -> barplot_relative('p__Cyanobacteria', 'o', 5, kraken_mpa_table)
    """
    
    fig, (ax_1, ax_2) = plt.subplots(2, 1, figsize=(16, 8), gridspec_kw={'height_ratios': [1, 0.3], 'hspace': 0.6})
    tax_level_table = taxa_table.loc[[str(i) for i in taxa_table.index if taxon_name in str(i) and i.split("|")[-1].startswith(f"{taxa_level}")]]
  
    top_n_by_sum = tax_level_table.loc[tax_level_table.sum(axis=1).nlargest(top_n).index]
    others_row = tax_level_table.loc[~tax_level_table.index.isin(top_n_by_sum.index)].sum()
    result_df = pd.concat([top_n_by_sum, pd.DataFrame(others_row).T.rename(index={0: 'Others'})])
    
    tax_level_table = result_df / result_df.sum()
    tax_level_table.index = ["|".join([j for j in i.split("|") if "x__" not in j]) for i in tax_level_table.index]
    
    tax_level_table.T.plot(kind='bar', stacked=True, legend=False, color=DISTINCT_COLORS, ax=ax_1, edgecolor='k', width=0.8)    
    ax_2.set_axis_off()
    ax_2.set_line_color = 'black'
    ax_1.set_ylabel('Relative abundance')
    ax_1.legend(
            loc='center left', 
            bbox_to_anchor=(1, 0.5), 
            title=''
        )
    fig.suptitle(f"Relative abundance of top {top_n} {taxon_name} at {taxa_level} level")
    fig.subplots_adjust(right=0.4)
    plt.tight_layout()
    plt.show()


def main():

    TAX_LEVELS_MAP = {'p': 6, 'c': 5, 'o': 4, 'f': 3, 'g': 2, 's': 1}
    PREFIXES = ["p__", "c__", "o__", "f__", "g__", "s__"]

    parser = argparse.ArgumentParser(description="A script to produce stacked barplot of microboal communities.")
    
    parser.add_argument("--input", type=str, required=True, help="Input Kraken MPA (MetaPhlAn format) file.")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--barplot_relative_to_all', action='store_true',
                       help="Barplot relative to all. Requires --tax-level and --top_n.")
    group.add_argument('--barplot_relative_to_parent', action='store_true',
                       help="Barplot relative to a parent. Requires --tax-name, --tax-level and --top_n.")
    
    parser.add_argument('--tax_name', type=str, help="Specify the exact name of the parent taxa (ex: p__Cyanobacteria).")
    parser.add_argument('--tax_level', type=str, help="Specify the taxonomic level [p/c/o/f/g/s].")
    parser.add_argument('--top_n', type=int, help="Specify the top N taxa to display in colors (max 20).")

    args = parser.parse_args()

    main_df = pd.read_csv(args.input, sep="\t", index_col=0)

    if args.tax_level:
        if args.tax_level not in TAX_LEVELS_MAP:
            parser.error("--tax_level must be one of the following: p, c, o, f, g, s.")

    if args.tax_name:
        if not any(args.tax_name.startswith(prefix) for prefix in PREFIXES):
            parser.error("--tax_name must start with one of the following prefixes: p__, c__, o__, f__, g__, s__.")

    if args.top_n:
        if args.top_n > 20:
            parser.error("--top_n must be max. 20.")

    if args.barplot_relative_to_all:
        if not args.tax_level or not args.top_n:
            parser.error("--barplot_relative_to_all requires --tax_level and --top_n.")
    
    if args.barplot_relative_to_parent:
        if not args.tax_name or not args.tax_level or not args.top_n:
            parser.error("--barplot_relative_to_parent requires --tax_name, --tax_level and --top_n.")

        if args.tax_name and args.tax_level:
            if TAX_LEVELS_MAP[args.tax_name.split("__")[0]] <= TAX_LEVELS_MAP[args.tax_level]:
                parser.error("--tax_level must be an inferior rank than --tax_name.")

    if args.barplot_relative_to_all:
        barplot_relative_to_all(args.tax_level, args.top_n, main_df)

    if args.barplot_relative_to_parent:
        barplot_relative_to_parent(args.tax_name, args.tax_level, args.top_n, main_df) 


if __name__ == "__main__":
    main()
