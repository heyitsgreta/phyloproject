# Phylogenetics Project Plan

## Context
University project for a bioinformatics/phylogenetics course. The goal is to analyze an unknown protein query sequence against a database of 539 TSR3 protein sequences from vertebrates using BLAST, MSA, and phylogenetic tree construction -- all in R. The project requires a well-documented R script with AI usage disclosed.

## Files
- **Input:** `context/database.fasta` (539 sequences, JSON metadata in headers), `context/query_sequence.fasta` (unknown ~300aa protein)
- **Output:** `phylo_analysis.R` (main script), `output/` directory (plots, BLAST results, alignment files)
- **Existing:** `hello.R` (can be ignored)

## Prerequisites
- Install BLAST: `brew install blast`
- R packages: Biostrings, seqinr, stringr, msa, ggplot2, ggmsa, pwalign, gplots, dendextend, phytools, ape, jsonlite

---

## Implementation Steps

### Step 0: Preamble
- Header comment documenting AI usage (required by assignment)
- Load all libraries, define paths, create `output/` and `blastdb/` directories
- Check that `blastp`/`makeblastdb` are available via `Sys.which()`

### Step 1: BLAST Search
1. **Build BLAST DB** from `database.fasta` using `system("makeblastdb ...")`
2. **Run BLASTP** with output format `"7 delim= sacc length pident gaps score bitscore evalue"`
3. **Parse results** -- filter `#` comment lines, read into data frame, sort by e-value
4. **Extract top 30 hits** by accession
5. **Evaluate substitution matrix** -- check median % identity of top 30:
   - `>80%` → BLOSUM80, `<40%` → BLOSUM45, else keep BLOSUM62
6. **Re-run BLAST** if matrix changed, re-extract top 30

### Step 2: Multiple Sequence Alignment
1. **Extract top 30 sequences** from `database.fasta` using `Biostrings::readAAStringSet()` and matching accessions
2. **Rename to binomial nomenclature** -- parse JSON from FASTA headers with `jsonlite::fromJSON()`, extract `organism_name`, take first two words as `Genus_species`, disambiguate duplicates with `_1`/`_2` suffix
3. **Add query sequence** as `Query_Unknown`
4. **Run MSA** with `msa(seqs, method = "ClustalW")`
5. **Visualize** with `ggmsa` (show a representative window ~100 positions) and save PDF

### Step 3: Distance Matrix & Heatmap
1. **Write alignment** to FASTA, read with `seqinr::read.alignment()`
2. **Compute distance matrix** with `dist.alignment(aln, matrix = "identity")`
3. **Heatmap** with `gplots::heatmap.2()` -- blue-white-red palette, save PDF
4. **Print summary** statistics (range, mean distance)

### Step 4: Tree Properties
1. **Four-point condition** -- for all C(31,4) = 31,465 quartets, check if the two largest of three pairwise sums are equal
2. **Ultrametric property** -- for all C(31,3) = 4,495 triples, check if the two largest distances are equal
3. Print violation counts and interpretation

### Step 5: Phylogenetic Trees
1. **UPGMA** via `hclust(dist, method = "average")` → `as.phylo()`
2. **Neighbor-Joining** via `ape::nj(dist)`
3. **Additional trees** -- complete linkage, Ward's method
4. **Visualize** all trees as PDFs
5. **Compare** -- tanglegram via `dendextend`, cophenetic correlation
6. **Bootstrap** (100 replicates) -- resample alignment columns, recompute distance + NJ tree, plot support values with `nodelabels()`
7. **Discussion comments** -- biological interpretation, comparison with vertebrate taxonomy, limitations

---

## Verification
1. Run `Rscript phylo_analysis.R` and confirm it completes without errors
2. Check `output/` directory contains: `blast_results.txt`, `alignment.fasta`, `msa_visualization.pdf`, `distance_heatmap.pdf`, `upgma_tree.pdf`, `nj_tree.pdf`, `tanglegram.pdf`, `nj_bootstrap.pdf`
3. Verify BLAST results are parsed correctly (30 rows in top hits)
4. Verify sequence names are in `Genus_species` format
