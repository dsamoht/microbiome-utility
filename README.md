# varied scripts for varied tasks
## [collapse_asv.py](https://github.com/dsamoht/utility/blob/main/collapse_asv.py)
- __dependencies__: pandas, [cdhit](https://github.com/weizhongli/cdhit)  
- __info__: This script aims to speed up DADA2's `collapseNoMismatch` by using `cd-hit-est` with 100% identity clustering. ASVs of different lengths but of 100% identity on the shorter sequence are collapsed (summed). The most abundant ASV across samples is kept as the representative.  
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
