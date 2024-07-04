# motifscanR

Scan input genomic regions with known DNA motifs

Given a set of input genomic regions, `motifscanR` scans the sequences to detect the occurrences of known motifs. It can also applies a statistical test on each motif to check whether the motif is significantly over- or under-represented (enriched or depleted) in the input genomic regions compared to another set of control regions.

# Citation

To cite the `motifscanR` package in publications, please use

> [Sun, H., Wang, J., Gong, Z. et al.  *Quantitative integration of epigenomic variation and transcription factor binding using MAmotif toolkit identifies an important role of IRF2 as transcription activator at gene promoters.* Cell Discov 4, **38** (2018).](https://doi.org/10.1038/s41421-018-0045-y)

# Installation

The latest version release of `motifscanR` could install with:


```r
BiocManager::install('LouisKwok-PICB/motifscanR')
```

# Usage

`motifscanR` could be used for genomic regions motif enrichment analysis
with `motifScan` function. 
```r
library(motifscanR)
library(BSgenome.Hsapiens.UCSC.hg19)

# Get a motif pwms
example_motifs <- getJasparMotifs(species = "Homo sapiens",
                                  collection = "CORE")

# Make a set of peaks
peaks <- GenomicRanges::GRanges(seqnames = c("chr1","chr2","chr2"),
                ranges = IRanges::IRanges(start = c(76585873,42772928,100183786),
                                          width = 500))

# Scan motif for example motifs
motif_ix <- motifScan(example_motifs, peaks, genome = "BSgenome.Hsapiens.UCSC.hg19")
```
The input object of genomic regions could be either [GenomicRanges](https://kasperdanielhansen.github.io/genbioconductor/html/GenomicRanges_GRanges.html), [DNAStringSet](https://kasperdanielhansen.github.io/genbioconductor/html/Biostrings.html), 
[DNAString](https://kasperdanielhansen.github.io/genbioconductor/html/Biostrings.html), or character vector. You could see more detail with `?motifScan`

`motifscanR` could also be used for genomic regions motif enrichment analysis
 between two sets of genomic regions with `motifEnrichment` function. 

```r
# Get a motif pwms
example_motifs <- getJasparMotifs(species = "Homo sapiens",
                                  collection = "CORE")
# Make a set of input regions
Input <- GenomicRanges::GRanges(seqnames = c("chr1","chr2","chr2"),
                                ranges = IRanges::IRanges(start = c(76585873,42772928,100183786),
                                                          width = 500))
# Make a set of control regions
Control <- GenomicRanges::GRanges(seqnames = c("chr1","chr3","chr5"),
                                  ranges = IRanges::IRanges(start = c(453123,6524593,100184233),
                                                            width = 500))
# Scan motif for example motifs
motif_ix_input <- motifScan(example_motifs, Input, genome = "BSgenome.Hsapiens.UCSC.hg19")
motif_ix_control <- motifScan(example_motifs, Control, genome = "BSgenome.Hsapiens.UCSC.hg19")

# Find Enrichment motif of input by control
Enrichment_result <- motifEnrichment(motif_ix_input, motif_ix_control)
```
You could type `?motifEnrichment` for more detail.
