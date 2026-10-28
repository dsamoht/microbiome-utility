## varied scripts for varied bioinformatics tasks in microbiome study, in no particular order
## [kreport_mpa_barplot.py](https://github.com/dsamoht/utility/blob/main/kreport_mpa_barplot.py)
- __info__: This script can produce 2 types of stacked barplots from Kraken classification reports and also plot dendrogram based on Bray-Curtis distance:  
  1. Top *n* taxa at a specific taxonomic level relative to all classified sequences (ex: barplot of top 10 phyla)
  2. Top *n* taxa at a specific taxonomic level relative to a parent taxa (ex: barplot of top 10 cyanobacterial genus inside the *Cyanobacteria* phylum)
- __input__: multiple kraken reports merged into a single table via 2 commands:  
  1. KrakenTools' kreport2mpa.py (with --intermediate-ranks parameters)
  2. KrakenTools' combine_mpa.py
- __dependencies__: pandas, matplotlib, scipy  
- __usage__:
```
python kreport_mpa_barplot.py -i [KREPORT_MPA] -l [p/c/o/f/g/s] -n [INT]
python kreport_mpa_barplot.py -i [KREPORT_MPA] -l [p/c/o/f/g/s] -n [INT] -p [PARENT_TAXA_STRING]
```
  - `-i/--input`: combined Kraken MPA table
  - `-l/--level`: target taxonomic level (p/c/o/f/g/s)
  - `-n/--top`: number of top taxa to display in colors (default 10, max 20)
  - `-p/--parent`: optional parent taxa (ex: `p__Cyanobacteriota`); abundances are then relative to that parent

  Plot top 5 orders within Cyanobacteria:
```
python kreport_mpa_barplot.py -i kraken_mpa.tsv -l o -n 5 -p p__Cyanobacteriota
```
- __example__: Top 15 prokaryotic phyla (input filtered to Bacteria and Archaea)
![alt text](img/barplot_example.png)



## [collapse_asv.py](https://github.com/dsamoht/utility/blob/main/collapse_asv.py)
- __dependencies__: pandas 
- __info__: This script aims to speed up DADA2's `collapseNoMismatch` by using python string matching. ASVs of different lengths but of 100% identity on the shorter sequence are collapsed (summed). The most abundant ASV across samples is kept as the representative.  
- __usage__:
```
python collapse_asv.py [SEQTAB] [COLLAPSED_SEQTAB]
```
- __input format__:  
[SEQTAB] must be comma-separated (ASV as column, sample as row, sample name as 1st column):  
(This format is typically the output of R write.csv() applied on a DADA2 seqtab object)
```
,CACGGA,ACACG,ATACCG,...    
sample_1,0,0,961,...    
sample_2,4,0,5,...    
sample_3,1,6,78,...    
...
```
## [fasta_sampler.py](https://github.com/dsamoht/utility/blob/main/fasta_sampler.py)
- __info__: This script samples a fasta file from a list of headers given in entry.
- __usage__:
```
python fasta_sampler.py [FASTA] [HEADERS] > [SAMPLED_FASTA]
```
- __input format__:  
[FASTA] : standard fasta  
[HEADERS] : file with 1 header of interest per line (with or without ">")
