###############################################################################
# Phylogenetic Analysis of TSR3 Protein Sequences
#
# Course: Bioinformatics / Phylogenetics
# Date: 2026-02-17
#
# AI Usage Disclosure:
#   This script was developed with assistance from Claude (Anthropic).
#   Claude was used to help structure the analysis pipeline, write R code for
#   BLAST integration, multiple sequence alignment, distance matrix computation,
#   tree property checks, and phylogenetic tree construction/comparison.
#   All code was reviewed and understood before submission.
#
# Description:
#   Analyzes an unknown protein query sequence against a database of 539 TSR3
#   protein sequences from vertebrates. Pipeline: BLAST search -> MSA ->
#   distance matrix -> tree property analysis -> phylogenetic tree construction.
###############################################################################

# ==== Step 0: Preamble ====

# Load libraries
library(Biostrings)
library(seqinr)
library(stringr)
library(msa)
library(ggplot2)
library(ggmsa)
library(pwalign)
library(gplots)
library(dendextend)
library(phytools)
library(ape)
library(jsonlite)

# Define paths
project_dir   <- getwd()
context_dir   <- file.path(project_dir, "context")
output_dir    <- file.path(project_dir, "output")
blastdb_dir   <- file.path(project_dir, "blastdb")
db_fasta      <- file.path(context_dir, "database.fasta")
query_fasta   <- file.path(context_dir, "query_sequence.fasta")

# Create directories
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(blastdb_dir, showWarnings = FALSE, recursive = TRUE)

# Check BLAST availability
blastp_path    <- Sys.which("blastp")
makeblastdb_path <- Sys.which("makeblastdb")
if (blastp_path == "" || makeblastdb_path == "") {
  stop("BLAST tools (blastp/makeblastdb) not found in PATH. Install with: brew install blast")
}
cat("BLAST tools found:\n  blastp:", blastp_path, "\n  makeblastdb:", makeblastdb_path, "\n\n")

# ==== Step 1: BLAST Search ====
cat("==== Step 1: BLAST Search ====\n")

# Build BLAST database
db_name <- file.path(blastdb_dir, "tsr3_db")
cmd_makedb <- sprintf(
  '%s -in "%s" -dbtype prot -out "%s" -parse_seqids',
  makeblastdb_path, db_fasta, db_name
)
cat("Building BLAST database...\n")
system(cmd_makedb, ignore.stdout = TRUE, ignore.stderr = TRUE)

# Function to run BLAST and parse results
run_blast <- function(query, db, matrix = "BLOSUM62", outfile = NULL) {
  if (is.null(outfile)) outfile <- file.path(output_dir, "blast_results.txt")

  cmd <- sprintf(
    '%s -query "%s" -db "%s" -matrix %s -outfmt "7 sacc length pident gaps score bitscore evalue" -out "%s"',
    blastp_path, query, db, matrix, outfile
  )
  system(cmd)

  # Parse results: skip comment lines starting with #
  lines <- readLines(outfile)
  data_lines <- lines[!grepl("^#", lines)]
  data_lines <- data_lines[nchar(trimws(data_lines)) > 0]

  if (length(data_lines) == 0) {
    stop("No BLAST hits found!")
  }

  df <- read.table(text = data_lines, sep = "\t", stringsAsFactors = FALSE,
                   col.names = c("sacc", "length", "pident", "gaps", "score", "bitscore", "evalue"))
  df$evalue <- as.numeric(df$evalue)
  df$pident <- as.numeric(df$pident)

  # Sort by e-value ascending
  df <- df[order(df$evalue), ]
  return(df)
}

# Initial BLAST run with default BLOSUM62
cat("Running initial BLASTP with BLOSUM62...\n")
blast_results <- run_blast(query_fasta, db_name, matrix = "BLOSUM62")
cat(sprintf("Total hits: %d\n", nrow(blast_results)))

# Extract top 30 unique accessions
top_acc <- unique(blast_results$sacc)[1:30]
top30 <- blast_results[blast_results$sacc %in% top_acc, ]
top30 <- top30[!duplicated(top30$sacc), ]
top30 <- top30[order(top30$evalue), ]

# Evaluate substitution matrix based on median % identity
median_pident <- median(top30$pident)
cat(sprintf("Median %% identity of top 30: %.1f%%\n", median_pident))

chosen_matrix <- "BLOSUM62"
if (median_pident > 80) {
  chosen_matrix <- "BLOSUM80"
} else if (median_pident < 40) {
  chosen_matrix <- "BLOSUM45"
}

# Re-run BLAST if matrix changed
if (chosen_matrix != "BLOSUM62") {
  cat(sprintf("Re-running BLASTP with %s...\n", chosen_matrix))
  blast_results <- run_blast(query_fasta, db_name, matrix = chosen_matrix)
  top_acc <- unique(blast_results$sacc)[1:30]
  top30 <- blast_results[blast_results$sacc %in% top_acc, ]
  top30 <- top30[!duplicated(top30$sacc), ]
  top30 <- top30[order(top30$evalue), ]
} else {
  cat("Keeping BLOSUM62 (median identity in 40-80% range).\n")
}

cat(sprintf("Final matrix: %s\n", chosen_matrix))
cat(sprintf("Top 30 hits (sorted by e-value):\n"))
print(top30[, c("sacc", "pident", "evalue", "bitscore")], row.names = FALSE)
cat("\n")

# ==== Step 2: Multiple Sequence Alignment ====
cat("==== Step 2: Multiple Sequence Alignment ====\n")

# Read all database sequences
all_seqs <- readAAStringSet(db_fasta)

# Match top 30 accessions to sequence names
# Accessions in BLAST output correspond to the part before the first space
seq_acc <- sub(" .*", "", names(all_seqs))
matched_idx <- match(top_acc, seq_acc)
matched_idx <- matched_idx[!is.na(matched_idx)]

if (length(matched_idx) < 30) {
  cat(sprintf("Warning: Only matched %d of 30 accessions\n", length(matched_idx)))
}

top_seqs <- all_seqs[matched_idx]

# Rename to binomial nomenclature from JSON metadata in headers
rename_to_binomial <- function(seqs) {
  full_names <- names(seqs)
  binomials <- character(length(full_names))

  for (i in seq_along(full_names)) {
    header <- full_names[i]
    # Extract JSON from header
    json_start <- regexpr("\\{", header)
    if (json_start > 0) {
      json_str <- substring(header, json_start)
      meta <- tryCatch(fromJSON(json_str), error = function(e) NULL)
      if (!is.null(meta) && !is.null(meta$organism_name)) {
        words <- strsplit(meta$organism_name, " ")[[1]]
        binomials[i] <- paste(words[1:min(2, length(words))], collapse = "_")
      } else {
        binomials[i] <- paste0("Unknown_", i)
      }
    } else {
      binomials[i] <- paste0("Unknown_", i)
    }
  }

  # Disambiguate duplicates by appending _1, _2, etc.
  dup_table <- table(binomials)
  dup_names <- names(dup_table[dup_table > 1])
  for (dn in dup_names) {
    idx <- which(binomials == dn)
    for (j in seq_along(idx)) {
      binomials[idx[j]] <- paste0(dn, "_", j)
    }
  }

  names(seqs) <- binomials
  return(seqs)
}

top_seqs <- rename_to_binomial(top_seqs)
cat("Top 30 species:\n")
cat(paste(" ", names(top_seqs)), sep = "\n")
cat("\n")

# Add query sequence
query_seq <- readAAStringSet(query_fasta)
names(query_seq) <- "Query_Unknown"

# Combine: query + top 30
all_for_msa <- c(query_seq, top_seqs)
cat(sprintf("Total sequences for MSA: %d\n", length(all_for_msa)))

# Run MSA with ClustalW
cat("Running ClustalW alignment...\n")
msa_result <- msa(all_for_msa, method = "ClustalW")
cat("MSA completed.\n")

# Visualize MSA with ggmsa
# Pick a representative window in the middle of the alignment
aln_width <- ncol(msa_result)
window_start <- max(1, round(aln_width / 2) - 50)
window_end <- min(aln_width, window_start + 99)

# Convert to AAStringSet for ggmsa
msa_ss <- AAStringSet(msa_result)

p_msa <- ggmsa(msa_ss, start = window_start, end = window_end,
               char_width = 0.5, seq_name = TRUE) +
  ggtitle(sprintf("MSA Visualization (positions %d-%d)", window_start, window_end))

ggsave(file.path(output_dir, "msa_visualization.pdf"), plot = p_msa,
       width = 16, height = 12)
cat("MSA visualization saved to output/msa_visualization.pdf\n\n")

# ==== Step 3: Distance Matrix & Heatmap ====
cat("==== Step 3: Distance Matrix & Heatmap ====\n")

# Write alignment to FASTA for seqinr
aln_file <- file.path(output_dir, "alignment.fasta")
writeXStringSet(msa_ss, filepath = aln_file)

# Read alignment with seqinr
aln <- read.alignment(file = aln_file, format = "fasta")

# Compute distance matrix (identity-based)
dist_mat <- dist.alignment(aln, matrix = "identity")
dist_matrix <- as.matrix(dist_mat)

cat(sprintf("Distance matrix dimensions: %d x %d\n", nrow(dist_matrix), ncol(dist_matrix)))
cat(sprintf("Distance range: %.4f - %.4f\n", min(dist_mat), max(dist_mat)))
cat(sprintf("Mean distance: %.4f\n", mean(dist_mat)))
cat("\n")

# Heatmap
pdf(file.path(output_dir, "distance_heatmap.pdf"), width = 14, height = 12)
heatmap.2(dist_matrix,
          col = colorRampPalette(c("blue", "white", "red"))(100),
          trace = "none",
          margins = c(12, 12),
          main = "Pairwise Distance Heatmap (Identity-Based)",
          key.title = "Distance",
          cexRow = 0.7, cexCol = 0.7,
          symm = TRUE)
dev.off()
cat("Heatmap saved to output/distance_heatmap.pdf\n\n")

# ==== Step 4: Tree Properties ====
cat("==== Step 4: Tree Properties ====\n")

n_taxa <- nrow(dist_matrix)
labels <- rownames(dist_matrix)

# Four-point condition check
# For an additive tree, for any 4 taxa, the two largest of the three pairwise
# sums d(i,j)+d(k,l), d(i,k)+d(j,l), d(i,l)+d(j,k) should be equal.
cat("Checking four-point condition...\n")
quartets <- combn(n_taxa, 4)
n_quartets <- ncol(quartets)
fp_violations <- 0
tolerance <- 1e-6

for (q in seq_len(n_quartets)) {
  idx <- quartets[, q]
  i <- idx[1]; j <- idx[2]; k <- idx[3]; l <- idx[4]

  s1 <- dist_matrix[i, j] + dist_matrix[k, l]
  s2 <- dist_matrix[i, k] + dist_matrix[j, l]
  s3 <- dist_matrix[i, l] + dist_matrix[j, k]

  sums <- sort(c(s1, s2, s3))
  # The two largest should be equal for additive tree
  if (abs(sums[3] - sums[2]) > tolerance) {
    fp_violations <- fp_violations + 1
  }
}
cat(sprintf("Four-point condition: %d / %d quartets violated (%.1f%%)\n",
            fp_violations, n_quartets, 100 * fp_violations / n_quartets))
if (fp_violations == 0) {
  cat("  -> Distance matrix is perfectly additive (tree metric).\n")
} else {
  cat("  -> Distance matrix is NOT perfectly additive.\n")
  cat("  -> This is expected for real biological data due to unequal evolutionary rates.\n")
}

# Ultrametric property check
# For an ultrametric tree, for any 3 taxa, the two largest distances should be equal.
cat("\nChecking ultrametric property...\n")
triples <- combn(n_taxa, 3)
n_triples <- ncol(triples)
um_violations <- 0

for (t in seq_len(n_triples)) {
  idx <- triples[, t]
  i <- idx[1]; j <- idx[2]; k <- idx[3]

  dists <- sort(c(dist_matrix[i, j], dist_matrix[i, k], dist_matrix[j, k]))
  # The two largest should be equal for ultrametric
  if (abs(dists[3] - dists[2]) > tolerance) {
    um_violations <- um_violations + 1
  }
}
cat(sprintf("Ultrametric property: %d / %d triples violated (%.1f%%)\n",
            um_violations, n_triples, 100 * um_violations / n_triples))
if (um_violations == 0) {
  cat("  -> Distance matrix is ultrametric (consistent with molecular clock).\n")
} else {
  cat("  -> Distance matrix is NOT ultrametric.\n")
  cat("  -> This suggests unequal rates of evolution among lineages (no strict molecular clock).\n")
}
cat("\n")

# ==== Step 5: Phylogenetic Trees ====
cat("==== Step 5: Phylogenetic Trees ====\n")

# UPGMA tree
cat("Building UPGMA tree...\n")
upgma_hclust <- hclust(dist_mat, method = "average")
upgma_tree <- as.phylo(upgma_hclust)

pdf(file.path(output_dir, "upgma_tree.pdf"), width = 12, height = 10)
plot(upgma_tree, main = "UPGMA Tree", cex = 0.7)
dev.off()
cat("UPGMA tree saved to output/upgma_tree.pdf\n")

# Neighbor-Joining tree
cat("Building Neighbor-Joining tree...\n")
nj_tree <- nj(dist_mat)

pdf(file.path(output_dir, "nj_tree.pdf"), width = 12, height = 10)
plot(nj_tree, main = "Neighbor-Joining Tree", cex = 0.7)
dev.off()
cat("NJ tree saved to output/nj_tree.pdf\n")

# Additional trees: Complete linkage and Ward's method
cat("Building additional trees (complete linkage, Ward)...\n")
complete_hclust <- hclust(dist_mat, method = "complete")
ward_hclust <- hclust(dist_mat, method = "ward.D2")

pdf(file.path(output_dir, "complete_linkage_tree.pdf"), width = 12, height = 10)
plot(as.phylo(complete_hclust), main = "Complete Linkage Tree", cex = 0.7)
dev.off()

pdf(file.path(output_dir, "ward_tree.pdf"), width = 12, height = 10)
plot(as.phylo(ward_hclust), main = "Ward's Method Tree", cex = 0.7)
dev.off()

# Tanglegram: UPGMA vs Complete Linkage
# (NJ trees are not ultrametric and cannot be converted to hclust/dendrogram directly)
cat("Creating tanglegram (UPGMA vs Complete Linkage)...\n")
upgma_dend <- as.dendrogram(upgma_hclust)
complete_dend <- as.dendrogram(complete_hclust)

# Untangle for better visualization
dend_list <- dendlist(upgma_dend, complete_dend)
dend_list <- dendextend::untangle(dend_list, method = "step2side")

pdf(file.path(output_dir, "tanglegram.pdf"), width = 14, height = 10)
tanglegram(dend_list,
           main_left = "UPGMA", main_right = "Complete Linkage",
           margin_inner = 7, lwd = 1.5, cex_main = 1.2,
           common_subtrees_color_branches = TRUE)
dev.off()
cat("Tanglegram saved to output/tanglegram.pdf\n")

# Also create a cophyloplot comparing UPGMA and NJ as phylo objects
pdf(file.path(output_dir, "cophylo_upgma_nj.pdf"), width = 14, height = 10)
obj <- cophylo(upgma_tree, nj_tree)
plot(obj, mar = c(0.5, 0.5, 4, 0.5), cex = 0.6)
title("UPGMA vs Neighbor-Joining Comparison")
dev.off()
cat("Cophyloplot (UPGMA vs NJ) saved to output/cophylo_upgma_nj.pdf\n")

# Cophenetic correlation
coph_upgma <- cophenetic(upgma_hclust)
coph_nj <- cophenetic(nj_tree)

# Align labels for correlation
common_labels <- intersect(rownames(as.matrix(coph_upgma)), rownames(as.matrix(coph_nj)))
coph_upgma_mat <- as.matrix(coph_upgma)[common_labels, common_labels]
coph_nj_mat <- as.matrix(coph_nj)[common_labels, common_labels]

# Cophenetic correlation with original distances
orig_mat <- dist_matrix[common_labels, common_labels]
cor_upgma <- cor(as.dist(orig_mat), as.dist(coph_upgma_mat))
cor_nj <- cor(as.dist(orig_mat), as.dist(coph_nj_mat))
cat(sprintf("Cophenetic correlation (UPGMA vs original): %.4f\n", cor_upgma))
cat(sprintf("Cophenetic correlation (NJ vs original): %.4f\n", cor_nj))

# Bootstrap analysis (NJ tree, 100 replicates)
cat("\nBootstrap analysis (100 replicates)...\n")

# Read alignment as matrix of characters
aln_matrix <- base::do.call(base::rbind, base::strsplit(as.character(aln$seq), ""))
n_sites <- ncol(aln_matrix)
n_seqs <- nrow(aln_matrix)
seq_names <- aln$nam

bootstrap_trees <- list()
n_boot <- 100

for (b in seq_len(n_boot)) {
  # Resample columns with replacement
  boot_cols <- sample(n_sites, n_sites, replace = TRUE)
  boot_aln_mat <- aln_matrix[, boot_cols]

  # Reconstruct alignment as seqinr format
  boot_seqs <- apply(boot_aln_mat, 1, paste, collapse = "")
  boot_aln <- list(
    nb = n_seqs,
    nam = seq_names,
    seq = boot_seqs,
    com = NA
  )
  class(boot_aln) <- "alignment"

  # Compute distance and build NJ tree
  boot_dist <- tryCatch(
    dist.alignment(boot_aln, matrix = "identity"),
    error = function(e) NULL
  )

  if (!is.null(boot_dist) && all(is.finite(as.matrix(boot_dist)))) {
    boot_tree <- tryCatch(nj(boot_dist), error = function(e) NULL)
    if (!is.null(boot_tree)) {
      bootstrap_trees[[length(bootstrap_trees) + 1]] <- boot_tree
    }
  }
}

cat(sprintf("Successfully computed %d / %d bootstrap trees.\n",
            length(bootstrap_trees), n_boot))

# Compute bootstrap support using prop.clades
if (length(bootstrap_trees) > 0) {
  boot_support <- prop.clades(nj_tree, bootstrap_trees)
  boot_pct <- round(100 * boot_support / length(bootstrap_trees))

  pdf(file.path(output_dir, "nj_bootstrap.pdf"), width = 14, height = 12)
  plot(nj_tree, main = "Neighbor-Joining Tree with Bootstrap Support (100 replicates)",
       cex = 0.7)
  # Only label internal nodes with support values
  nodelabels(boot_pct, cex = 0.6, bg = "lightyellow", frame = "rect")
  dev.off()
  cat("Bootstrap tree saved to output/nj_bootstrap.pdf\n")
}

# ==== Discussion ====
cat("\n==== Discussion ====\n")
cat("
Biological Interpretation:
- The query sequence groups with TSR3 (ribosome biogenesis protein) homologs,
  a conserved protein involved in 18S rRNA modification in eukaryotes.
- The phylogenetic trees should broadly reflect vertebrate taxonomy, with
  major clades (mammals, birds, reptiles, fish) forming distinct groups.

Tree Comparison:
- UPGMA assumes a molecular clock (equal rates of evolution), producing an
  ultrametric tree. This assumption is often violated in real data.
- Neighbor-Joining does not assume a molecular clock and generally fits
  non-ultrametric distances better, as reflected by cophenetic correlations.
- The tanglegram visualizes topological differences between the two methods.

Limitations:
- Using only 30 taxa from 539 available limits resolution.
- Single-gene phylogenies may not reflect species phylogeny due to
  incomplete lineage sorting or horizontal gene transfer.
- ClustalW is a heuristic aligner; results may vary with different methods.
- Bootstrap support helps quantify confidence in tree topology, but 100
  replicates is a minimum; 1000+ would be more robust.
")

cat("\n==== Analysis Complete ====\n")
cat("Output files:\n")
cat(paste(" ", list.files(output_dir, full.names = FALSE)), sep = "\n")
cat("\n")
