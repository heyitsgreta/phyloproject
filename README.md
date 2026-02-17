# Phylogenetic Analysis of TSR3 Protein Sequences

University project for a bioinformatics/phylogenetics course. Analyzes an unknown protein query sequence against a database of 539 TSR3 (ribosome biogenesis protein) sequences from vertebrates.

## What it does

1. **BLAST search** -- identifies the top 30 most similar sequences from the database, with automatic substitution matrix selection based on sequence identity
2. **Multiple Sequence Alignment** -- aligns the top hits + query using ClustalW
3. **Distance matrix & heatmap** -- computes pairwise identity-based distances
4. **Tree property checks** -- tests the four-point condition (additivity) and ultrametric property
5. **Phylogenetic trees** -- builds UPGMA, Neighbor-Joining, complete linkage, and Ward's method trees with bootstrap support (100 replicates)

## Getting started

### Prerequisites

- **R** (>= 4.5)
- **BLAST+** command-line tools

Install BLAST on macOS:
```bash
brew install blast
```

### Install R packages

Open R and run:
```r
install.packages("BiocManager", repos = "https://cloud.r-project.org")
BiocManager::install(c("Biostrings", "msa", "pwalign", "ggmsa"), ask = FALSE)
install.packages(c("seqinr", "stringr", "ggplot2", "gplots", "dendextend",
                    "phytools", "ape", "jsonlite"),
                 repos = "https://cloud.r-project.org")
```

### Run the analysis

```bash
Rscript phylo_analysis.R
```

Results are written to the `output/` directory.

## Project structure

```
context/
  database.fasta        # 539 TSR3 protein sequences (vertebrates)
  query_sequence.fasta  # Unknown ~300aa query protein
phylo_analysis.R        # Main analysis script
output/
  blast_results.txt     # BLAST output
  alignment.fasta       # MSA in FASTA format
  msa_visualization.pdf # Alignment visualization
  distance_heatmap.pdf  # Pairwise distance heatmap
  upgma_tree.pdf        # UPGMA tree
  nj_tree.pdf           # Neighbor-Joining tree
  nj_bootstrap.pdf      # NJ tree with bootstrap support
  tanglegram.pdf        # UPGMA vs Complete Linkage comparison
  cophylo_upgma_nj.pdf  # UPGMA vs NJ comparison
```

## AI disclosure

This project was developed with assistance from Claude (Anthropic). See the header of `phylo_analysis.R` for details.
