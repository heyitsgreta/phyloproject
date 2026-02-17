# Code Walkthrough: `phylo_analysis.R`

This document is a detailed walkthrough of the `phylo_analysis.R` script for readers who are new to R and/or bioinformatics. It explains not just *what* the code does, but *why* each step exists and what biological or computational problem it solves.

---

## 1. Overview

The script implements a **phylogenetic analysis pipeline** -- a sequence of steps that starts with a mystery protein sequence and ends with evolutionary trees showing how that protein relates to known sequences.

**What goes in:**

- `context/query_sequence.fasta` -- a single unknown protein sequence (the "query")
- `context/database.fasta` -- 539 known TSR3 protein sequences from vertebrates

**What comes out (in `output/`):**

- BLAST search results (which database sequences are most similar to the query)
- A multiple sequence alignment (MSA) visualization
- A pairwise distance heatmap
- Several phylogenetic trees (UPGMA, Neighbor-Joining, Complete Linkage, Ward's method)
- A tanglegram and cophyloplot comparing tree methods
- A bootstrap-annotated tree showing confidence in the branching pattern

**The pipeline in plain English:**

1. Find the 30 most similar sequences to the query (BLAST)
2. Line them all up so corresponding amino acids are in the same columns (MSA)
3. Measure how different each pair of sequences is (distance matrix)
4. Check mathematical properties of those distances (tree property tests)
5. Build evolutionary trees from the distances and assess confidence (trees + bootstrap)

---

## 2. R Basics for Newcomers

Here are the R syntax patterns you will encounter in the script:

| Syntax | Meaning | Example |
|--------|---------|---------|
| `x <- value` | Assignment. Stores `value` in variable `x`. R's equivalent of `=` in most languages. | `project_dir <- getwd()` |
| `library(pkg)` | Load a package so its functions are available. Like `import` in Python. | `library(ape)` |
| `df$column` | Access a column in a data frame (table) by name. | `blast_results$evalue` |
| `df[rows, cols]` | Subset a data frame. `df[1:5, ]` = first 5 rows, all columns. | `top30[, c("sacc", "pident")]` |
| `function(args) { body }` | Define a function. | `run_blast <- function(query, db) { ... }` |
| `sprintf("text %s", var)` | Format a string, inserting variables. `%s` = string, `%d` = integer, `%.1f` = 1-decimal float. | `sprintf("Total hits: %d", nrow(df))` |
| `file.path(a, b)` | Join path segments with the OS-appropriate separator. Safer than pasting `/` manually. | `file.path(project_dir, "output")` |
| `c(x, y, z)` | Combine values into a vector (a 1-D list of the same type). | `c("sacc", "pident", "evalue")` |
| `for (x in seq) { }` | Loop over a sequence. | `for (i in seq_along(names)) { ... }` |
| `x %>% f()` | Pipe operator (not used here, but common in R). Passes `x` as the first argument to `f()`. | -- |
| `TRUE` / `FALSE` | Logical values (booleans). | `showWarnings = FALSE` |
| `NULL` | R's "nothing" value. Used for optional defaults and to signal absence. | `if (is.null(outfile))` |
| `base::strsplit()` | The `::` operator calls a function from a specific package. Used here to avoid conflicts when multiple packages define the same function name. | `base::do.call(...)` |

---

## 3. Packages: What and Why

The script loads 12 packages. Each solves a specific problem:

### Sequence I/O and Manipulation

| Package | What it does | Why the script needs it |
|---------|-------------|----------------------|
| **Biostrings** | Read, write, and manipulate biological sequences (DNA, RNA, protein). Part of Bioconductor. | Reading FASTA files (`readAAStringSet`), storing aligned sequences as `AAStringSet` objects, writing alignments back out (`writeXStringSet`). |
| **seqinr** | General-purpose sequence analysis toolkit. | Reading alignments in a format compatible with its `dist.alignment()` function, which computes pairwise sequence distances. |
| **pwalign** | Pairwise sequence alignment tools. Part of Bioconductor. | Provides alignment scoring infrastructure that other Bioconductor packages depend on (loaded for compatibility). |

### Multiple Sequence Alignment

| Package | What it does | Why the script needs it |
|---------|-------------|----------------------|
| **msa** | Runs multiple sequence alignment algorithms (ClustalW, ClustalOmega, MUSCLE) from within R. | The `msa()` function performs the ClustalW alignment of all 31 sequences (query + top 30 hits). |

### Plotting and Visualization

| Package | What it does | Why the script needs it |
|---------|-------------|----------------------|
| **ggplot2** | The most widely-used R plotting library. Builds plots by layering graphical components. | Provides the plotting foundation that `ggmsa` builds on top of. The `ggtitle()` and `ggsave()` functions come from here. |
| **ggmsa** | Visualizes multiple sequence alignments as publication-quality plots, with amino acids colored by chemical property. | Generates the MSA visualization showing aligned positions with color-coded residues. |
| **gplots** | Extended plotting functions including heatmaps. | The `heatmap.2()` function draws the distance matrix heatmap with dendrograms, color keys, and customization options. |

### Tree Building and Comparison

| Package | What it does | Why the script needs it |
|---------|-------------|----------------------|
| **ape** | The core R package for phylogenetics ("Analysis of Phylogenetics and Evolution"). | Builds Neighbor-Joining trees (`nj()`), converts hierarchical clusters to `phylo` objects (`as.phylo()`), computes bootstrap support (`prop.clades()`), and plots trees. |
| **phytools** | Extended phylogenetic tools, building on `ape`. | Provides `cophylo()` for creating face-to-face tree comparisons (cophyloplots) between UPGMA and NJ trees. |
| **dendextend** | Manipulate and compare tree-like structures (dendrograms). | Creates tanglegrams (`tanglegram()`) and untangles crossing lines (`untangle()`) for cleaner side-by-side tree comparisons. |

### Utilities

| Package | What it does | Why the script needs it |
|---------|-------------|----------------------|
| **stringr** | Consistent, readable string manipulation functions. | Available for string operations throughout the script. Provides cleaner alternatives to base R string functions. |
| **jsonlite** | Parse and generate JSON data. | The FASTA headers contain JSON metadata with organism names. `fromJSON()` extracts the organism name so sequences can be labeled with species names instead of cryptic accession IDs. |

---

## 4. Section-by-Section Walkthrough

### Step 0: Preamble (Lines 20-54)

**What it does:** Loads packages, defines file paths, creates output directories, and checks that BLAST is installed.

**Why setup matters:**

R scripts don't automatically have access to every function -- you must explicitly load each package with `library()`. This is different from languages where a standard library is always available. Each `library()` call makes that package's functions available in the current session.

**How `dir.create()` works:**

```r
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
```

- `showWarnings = FALSE` -- don't print a warning if the directory already exists
- `recursive = TRUE` -- create parent directories if needed (like `mkdir -p` in the terminal)

**What `Sys.which()` does:**

```r
blastp_path <- Sys.which("blastp")
```

`Sys.which()` searches your system's PATH for an executable, the same thing that happens when you type a command in the terminal. It returns the full path (e.g., `/usr/local/bin/blastp`) or an empty string `""` if the tool isn't installed.

This check is important because the script depends on external BLAST tools that are *not* R packages -- they are separate command-line programs. If they are missing, the script stops immediately with a clear error message via `stop()`, rather than failing cryptically later.

---

### Step 1: BLAST Search (Lines 56-134)

**The biological problem:**

You have a mystery protein sequence and a database of 539 known TSR3 sequences. You need to find which database sequences are most similar to the query. This is a classic "needle in a haystack" problem, and BLAST (Basic Local Alignment Search Tool) is the standard solution.

**Why BLAST instead of aligning everything?**

Aligning the query against all 539 sequences simultaneously would be computationally expensive and unnecessary. BLAST uses a clever shortcut: it breaks the query into short "words," finds database sequences that share those words, and only does a full alignment for promising matches. This makes it orders of magnitude faster.

**What a substitution matrix is and why it matters:**

```r
chosen_matrix <- "BLOSUM62"
if (median_pident > 80) {
  chosen_matrix <- "BLOSUM80"
} else if (median_pident < 40) {
  chosen_matrix <- "BLOSUM45"
}
```

A substitution matrix tells BLAST how to score amino acid matches and mismatches. Not all mismatches are equally bad -- replacing one amino acid with a chemically similar one (e.g., leucine with isoleucine) is less damaging than replacing it with a very different one (e.g., leucine with aspartic acid).

The BLOSUM family of matrices is derived from observed substitution patterns in real proteins:

- **BLOSUM80** -- tuned for highly similar sequences (>80% identity). More strict about mismatches.
- **BLOSUM62** -- the general-purpose default. Good for moderate similarity (40-80%).
- **BLOSUM45** -- tuned for divergent sequences (<40% identity). More forgiving of mismatches.

The script does a first pass with BLOSUM62, checks the median percent identity of the top hits, and re-runs with a more appropriate matrix if needed. This adaptive approach ensures the scoring reflects the actual evolutionary distance in the data.

**How R talks to command-line tools:**

```r
system(cmd_makedb, ignore.stdout = TRUE, ignore.stderr = TRUE)
```

`system()` runs a shell command from within R, exactly as if you typed it in the terminal. The script builds the command string with `sprintf()` (string formatting) and passes it to `system()`. The `ignore.stdout` and `ignore.stderr` arguments suppress terminal output to keep the R console clean.

**How data frames work for tabular results:**

```r
df <- read.table(text = data_lines, sep = "\t", stringsAsFactors = FALSE,
                 col.names = c("sacc", "length", "pident", "gaps", "score", "bitscore", "evalue"))
```

A data frame is R's equivalent of a spreadsheet or SQL table. Each column has a name and a type. `read.table()` parses the tab-separated BLAST output into this structured format so you can access columns by name (e.g., `df$evalue`), filter rows, sort, and perform calculations.

- `text = data_lines` -- read from an in-memory string instead of a file
- `sep = "\t"` -- columns are separated by tabs
- `stringsAsFactors = FALSE` -- treat text as text, not as categorical variables (an old R gotcha)
- `col.names = ...` -- give human-readable names to each column

**The BLAST output columns:**

| Column | Meaning |
|--------|---------|
| `sacc` | Subject accession -- the ID of the matching database sequence |
| `length` | Alignment length |
| `pident` | Percent identity -- what fraction of aligned positions are identical |
| `gaps` | Number of gap positions in the alignment |
| `score` | Raw alignment score |
| `bitscore` | Normalized score (comparable across searches) |
| `evalue` | Expected value -- how many hits this good you'd expect by chance. Lower = more significant. |

---

### Step 2: Multiple Sequence Alignment (Lines 136-224)

**Why alignment is needed before building a tree:**

Phylogenetic trees are built by comparing sequences position-by-position. But raw sequences have different lengths and may have insertions or deletions. An MSA inserts gap characters (`-`) so that homologous positions (positions that descended from the same ancestral position) line up in the same column.

Without alignment, you'd be comparing position 50 of one sequence with position 50 of another, even though they might correspond to completely different parts of the protein. Alignment ensures you're comparing like with like.

**What ClustalW does conceptually:**

```r
msa_result <- msa(all_for_msa, method = "ClustalW")
```

ClustalW is a progressive alignment algorithm:

1. It first compares all pairs of sequences to estimate how similar they are
2. It builds a rough "guide tree" from those similarities
3. It aligns the two most similar sequences first
4. It progressively adds more sequences to the alignment, following the guide tree from most to least similar

This progressive approach is a practical compromise -- aligning all sequences simultaneously would be computationally intractable, but aligning them in order of similarity gives good results.

**The JSON parsing trick:**

```r
json_start <- regexpr("\\{", header)
if (json_start > 0) {
  json_str <- substring(header, json_start)
  meta <- tryCatch(fromJSON(json_str), error = function(e) NULL)
```

The FASTA headers in the database contain metadata in JSON format (embedded in the header line after the accession). The function finds the opening `{`, extracts everything from there to the end, and parses it with `fromJSON()`. This is how the script gets the organism name for each sequence.

`tryCatch()` is R's try/catch mechanism -- if the JSON is malformed, it returns `NULL` instead of crashing the script.

**Why renaming matters:**

```r
names(seqs) <- binomials
```

FASTA accession IDs like `XP_123456789.1` are meaningless to a human reader. The renaming function extracts binomial nomenclature (e.g., `Homo_sapiens`) from the JSON metadata so the trees and plots have readable, biologically meaningful labels. Duplicate species get suffixes (`_1`, `_2`) to stay unique.

---

### Step 3: Distance Matrix & Heatmap (Lines 226-256)

**What a distance matrix represents:**

```r
dist_mat <- dist.alignment(aln, matrix = "identity")
```

A distance matrix is a table where each cell `(i, j)` contains a number representing how different sequence `i` is from sequence `j`. The `identity` method computes `1 - (fraction of identical positions)`, so:

- A distance of **0** means the sequences are identical
- A distance of **1** means no positions match

This gives you a single number summarizing the evolutionary divergence between any two sequences, which is what tree-building algorithms need as input.

**Why identity-based distance?**

There are more sophisticated evolutionary distance models (e.g., Jukes-Cantor, Kimura) that correct for multiple substitutions at the same position. Identity-based distance is simpler but works well for protein sequences with moderate divergence.

**What the heatmap shows biologically:**

```r
heatmap.2(dist_matrix, ...)
```

The heatmap is a visual summary of all pairwise distances. Sequences are reordered by hierarchical clustering (shown as dendrograms on the margins), so closely related sequences end up next to each other. Color encodes distance:

- **Blue** blocks = low distance = high similarity = recent common ancestor
- **Red** blocks = high distance = low similarity = distant relationship

You can instantly spot groups of closely related species (blue clusters along the diagonal) and identify which groups are most distant from each other.

---

### Step 4: Tree Properties (Lines 258-321)

**Why test tree properties before building trees?**

These tests check whether the distance data has mathematical properties that certain tree-building methods assume. If the data violates these assumptions, it tells you something biologically interesting and guides your choice of method.

**The four-point condition (additive tree test):**

```r
s1 <- dist_matrix[i, j] + dist_matrix[k, l]
s2 <- dist_matrix[i, k] + dist_matrix[j, l]
s3 <- dist_matrix[i, l] + dist_matrix[j, k]
sums <- sort(c(s1, s2, s3))
# The two largest should be equal for an additive tree
```

For any four sequences (a quartet), you can compute three sums of paired distances. If the distance data comes from a true tree, the two largest sums will always be equal. This is a mathematical property of tree-like (additive) distances.

**What it means biologically:** If all quartets pass, the distances fit a tree perfectly -- evolution along branches explains all the observed differences. Violations mean the data has "noise" that no single tree can perfectly explain. This is normal for real biological data due to factors like:

- Unequal rates of evolution in different lineages
- Multiple substitutions at the same site
- Measurement imprecision

**The ultrametric property:**

```r
dists <- sort(c(dist_matrix[i, j], dist_matrix[i, k], dist_matrix[j, k]))
# The two largest should be equal for ultrametric
```

For any three sequences, if the two largest pairwise distances are equal, the distances are **ultrametric**. Ultrametric distances arise when there is a strict molecular clock -- all lineages evolve at the same rate.

**What it means biologically:** If the data is ultrametric, a UPGMA tree (which assumes a molecular clock) is appropriate. If not (which is the usual case), you know that evolutionary rates vary across lineages, and methods like Neighbor-Joining that don't assume a clock will be more accurate.

---

### Step 5: Phylogenetic Trees (Lines 323-486)

**UPGMA (Unweighted Pair Group Method with Arithmetic Mean):**

```r
upgma_hclust <- hclust(dist_mat, method = "average")
```

UPGMA is a hierarchical clustering algorithm:

1. Start with each sequence as its own cluster
2. Find the two closest clusters and merge them
3. Recalculate distances to the new merged cluster (using the average of all pairwise distances)
4. Repeat until everything is in one cluster

UPGMA produces an **ultrametric** tree where the root-to-tip distance is the same for every sequence. This implicitly assumes a molecular clock. If evolutionary rates vary (which they usually do), UPGMA can produce misleading trees.

**Neighbor-Joining (NJ):**

```r
nj_tree <- nj(dist_mat)
```

Neighbor-Joining is a more flexible algorithm:

1. Compute a "corrected" distance for each sequence that accounts for how far it is from everything else
2. Find the pair of sequences whose joining minimizes the total tree length
3. Join them, creating a new internal node
4. Recalculate distances and repeat

NJ does **not** assume a molecular clock. It can produce trees with unequal branch lengths, which better represents reality when different lineages evolve at different rates. This is why NJ generally produces more accurate trees for real biological data.

**Why multiple methods?**

The script builds four trees (UPGMA, NJ, Complete Linkage, Ward's method) because no single method is guaranteed to recover the true evolutionary history. Comparing trees built with different algorithms reveals which relationships are robust (appear in all trees) and which are method-dependent (and thus less certain).

**What a tanglegram shows:**

```r
tanglegram(dend_list, main_left = "UPGMA", main_right = "Complete Linkage", ...)
```

A tanglegram draws two trees side-by-side with lines connecting the same species across both trees. If the trees agree, the lines run roughly parallel. If they disagree, the lines cross. The `untangle()` function reorders leaves to minimize crossings, making it easier to spot differences.

```r
dend_list <- dendextend::untangle(dend_list, method = "step2side")
```

Note the explicit `dendextend::untangle()` -- this uses the `::` syntax to call `untangle` specifically from the `dendextend` package, avoiding a conflict with another package that also defines an `untangle` function (see Section 5 on namespace conflicts).

**What cophenetic correlation measures:**

```r
cor_upgma <- cor(as.dist(orig_mat), as.dist(coph_upgma_mat))
cor_nj <- cor(as.dist(orig_mat), as.dist(coph_nj_mat))
```

The cophenetic distance between two sequences in a tree is the height of their lowest common ancestor (where their branches meet). Cophenetic correlation compares these tree-implied distances to the original distance matrix:

- A correlation of **1.0** means the tree perfectly represents the original distances
- Lower values mean the tree distorts some distances

This gives you a quantitative way to ask: "Which tree-building method did a better job of representing the actual sequence differences?" NJ typically wins for non-ultrametric data.

**What bootstrap does and why confidence matters:**

```r
boot_cols <- sample(n_sites, n_sites, replace = TRUE)
boot_aln_mat <- aln_matrix[, boot_cols]
```

Bootstrapping answers the question: "How confident are we in each branch of the tree?"

The procedure:

1. Take the original alignment (a matrix of sequences x positions)
2. Create a new "pseudo-alignment" by randomly sampling columns *with replacement* (some columns appear multiple times, others not at all)
3. Build a new tree from the resampled alignment
4. Repeat 100 times (or more)
5. For each branch in the original tree, count how many bootstrap trees contain the same branch

A bootstrap value of **95%** means that branch appeared in 95 out of 100 resampled trees -- it's strongly supported by the data. A value of **50%** means the branch is essentially a coin flip.

The key insight: if a branch is supported by many independent positions in the alignment, it will survive resampling. If it depends on just a few positions, different resamplings will produce different trees, and the bootstrap value will be low.

```r
boot_support <- prop.clades(nj_tree, bootstrap_trees)
```

`prop.clades()` from the `ape` package does the counting: for each internal node (branch point) in the original NJ tree, it counts how many bootstrap trees contain a branch that splits the same set of species.

---

## 5. Gotchas & Namespace Conflicts

When you load 12 packages, some of them define functions with the same name. R resolves this by letting the most recently loaded package "win" -- its version of the function masks the earlier one. This causes subtle bugs when you expect one version but get another.

The script handles three specific conflicts:

### `base::strsplit()` and `base::do.call()`

```r
aln_matrix <- base::do.call(base::rbind, base::strsplit(as.character(aln$seq), ""))
```

The `base::` prefix explicitly calls R's built-in versions of `strsplit()`, `do.call()`, and `rbind()`. Without it, a loaded package could have overridden these functions with its own versions that behave differently. Using the `base::` prefix guarantees you get the standard R behavior regardless of which packages are loaded.

### `dendextend::untangle()`

```r
dend_list <- dendextend::untangle(dend_list, method = "step2side")
```

The `phytools` package also exports a function called `untangle()`. Since `phytools` and `dendextend` are both loaded, whichever was loaded last would mask the other's `untangle()`. The `dendextend::` prefix ensures the correct version (the one that works on dendrograms) is called.

### General rule

If you see `package::function()` syntax in R code, it means the author is being defensive about namespace conflicts. This is good practice whenever you load many packages, especially from the Bioconductor ecosystem where function name collisions are common.

---

## Quick Reference: Output Files

| File | Produced by | Shows |
|------|------------|-------|
| `blast_results.txt` | Step 1 | Raw BLAST hits in tabular format |
| `msa_visualization.pdf` | Step 2 | Color-coded alignment of a representative window |
| `alignment.fasta` | Step 2 | Full alignment in FASTA format (used internally) |
| `distance_heatmap.pdf` | Step 3 | Clustered heatmap of all pairwise distances |
| `upgma_tree.pdf` | Step 5 | UPGMA phylogenetic tree |
| `nj_tree.pdf` | Step 5 | Neighbor-Joining phylogenetic tree |
| `complete_linkage_tree.pdf` | Step 5 | Complete Linkage tree |
| `ward_tree.pdf` | Step 5 | Ward's method tree |
| `tanglegram.pdf` | Step 5 | Side-by-side comparison of UPGMA vs Complete Linkage |
| `cophylo_upgma_nj.pdf` | Step 5 | Face-to-face comparison of UPGMA vs NJ |
| `nj_bootstrap.pdf` | Step 5 | NJ tree annotated with bootstrap support values |
