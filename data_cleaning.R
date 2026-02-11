################################################################################
#
#  Somatic Bridges and Sentinel Short-Forms:
#  A Network Analysis and In Silico Intervention Simulation of Internalizing
#  Comorbidity in 50,000 University Students (2019–2024)
#
#  COMPLETE ANALYSIS PIPELINE — H1 through H5
#  Lead Computational Social Scientist Script
#
#  Input:  scl90_data.csv
#  Output: Console report + saved figures and tables
#
################################################################################

# ──────────────────────────────────────────────────────────────────────────────
# 0.  ENVIRONMENT SETUP
# ──────────────────────────────────────────────────────────────────────────────

required_packages <- c(

  "bootnet",               # Network estimation (EBICglasso)
  "qgraph",                # Network visualization
  "igraph",                # Community detection (Walktrap, Louvain)
  "networktools",          # Bridge centrality, goldbricker
  "NetworkComparisonTest", # NCT for temporal invariance (H5)
  "dplyr",                 # Data wrangling
  "tidyr",                 # Data reshaping
  "tibble",                # Tibbles
  "purrr",                 # Functional programming
  "ggplot2",               # Plotting
  "pROC",                  # ROC analysis (H4)
  "effectsize",            # Cohen's d
  "mclust",                # Adjusted Rand Index
  "psych",                 # Descriptive statistics
  "readr",                 # CSV import
  "glue",                  # String interpolation
  "knitr"                  # Table formatting
)

install_if_missing <- function(pkg) {

  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
  }
}
invisible(lapply(required_packages, install_if_missing))
invisible(lapply(required_packages, library, character.only = TRUE))

cat("\n",
    strrep("=", 78), "\n",
    "  ANALYSIS PIPELINE INITIALISED — ",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
    strrep("=", 78), "\n\n")

set.seed(2024)  # Reproducibility

# ──────────────────────────────────────────────────────────────────────────────
# 1.  DATA IMPORT & ITEM DEFINITIONS
# ──────────────────────────────────────────────────────────────────────────────

cat("──── STAGE 0: Data Import ────\n")

raw <- read_csv("scl90_data.csv", show_col_types = FALSE)
cat(glue("  Raw sample loaded: N = {nrow(raw)}, variables = {ncol(raw)}"), "\n\n")

# --- Define the 35 internalizing items (SCL-90-R numbering) -----------------
# Somatization (12 items)
som_items <- paste0("Q", c(1, 4, 12, 27, 40, 42, 48, 49, 52, 53, 56, 58))
# Depression (13 items)
dep_items <- paste0("Q", c(5, 14, 15, 20, 22, 26, 29, 30, 31, 32, 54, 71, 79))
# Anxiety (10 items)
anx_items <- paste0("Q", c(2, 17, 23, 33, 39, 57, 72, 78, 80, 86))

all_items <- c(som_items, dep_items, anx_items)
stopifnot(length(all_items) == 35)

# Human-readable labels for reporting
item_labels <- c(
  Q1  = "Headaches",        Q4  = "Faintness",
  Q12 = "Chest pain",       Q27 = "Back pain",
  Q40 = "Nausea",           Q42 = "Muscle soreness",
  Q48 = "Temp dysreg",      Q49 = "Numbness",
  Q52 = "Throat lump",      Q53 = "Heavy limbs",
  Q56 = "Body weakness",    Q58 = "Death thoughts (Som)",
  Q5  = "Sexual disinterest", Q14 = "Low energy",
  Q15 = "Suicidal ideation",  Q20 = "Crying",
  Q22 = "Feeling trapped",    Q26 = "Self-blame",
  Q29 = "Loneliness",         Q30 = "Blue mood",
  Q31 = "Excessive worry",    Q32 = "Anhedonia",
  Q54 = "Hopelessness",       Q71 = "Death preoccup",
  Q79 = "Worthlessness",
  Q2  = "Nervousness",      Q17 = "Trembling",
  Q23 = "Sudden fear",      Q33 = "Fearfulness",
  Q39 = "Palpitations",     Q57 = "Tension",
  Q72 = "Panic",            Q78 = "Restlessness",
  Q80 = "Dread",            Q86 = "Fright cognitions"
)

# Subscale membership vector (for ARI comparisons)
subscale_membership <- setNames(

  c(rep("Somatization", 12), rep("Depression", 13), rep("Anxiety", 10)),
  all_items
)

# Binary community for bridge analysis: Somatic vs Affective (Dep + Anx)
bridge_communities <- setNames(
  c(rep("Somatic", 12), rep("Affective", 23)),
  all_items
)

# ──────────────────────────────────────────────────────────────────────────────
# 2.  DATA QUALITY ASSURANCE — THREE-STAGE SEQUENTIAL FILTERING
# ──────────────────────────────────────────────────────────────────────────────

cat("──── STAGE 1: Longstring Analysis ────\n")

# Check that all 35 items exist in the data
missing_items <- setdiff(all_items, names(raw))
if (length(missing_items) > 0) {
  stop(glue("Missing items in data: {paste(missing_items, collapse = ', ')}"))
}

# Ensure a cohort/year variable exists for H5
# Adapt column name as needed: expects 'cohort_year' or 'year'
if (!"cohort_year" %in% names(raw)) {
  if ("year" %in% names(raw)) {
    raw <- raw %>% rename(cohort_year = year)
  } else {
    warning("No cohort/year variable found. H5 temporal analysis will be skipped.")
    raw$cohort_year <- NA
  }
}

# Compute longstring: max consecutive identical responses across 35 items
compute_longstring <- function(row_vec) {
  r <- rle(row_vec)
  max(r$lengths)
}

item_matrix <- as.matrix(raw[, all_items])
raw$longstring <- apply(item_matrix, 1, compute_longstring)

n_before <- nrow(raw)
dat_s1 <- raw %>% filter(longstring < 15)
n_ls_excl <- n_before - nrow(dat_s1)
cat(glue("  Longstring ≥15 excluded: n = {n_ls_excl} ({round(n_ls_excl/n_before*100,1)}%)"), "\n")
cat(glue("  Remaining: N = {nrow(dat_s1)}"), "\n\n")

# --- Stage 2: Mahalanobis Distance Outlier Detection ------------------------
cat("──── STAGE 2: Mahalanobis Distance Outlier Detection ────\n")

item_mat_s1 <- as.matrix(dat_s1[, all_items])
mu <- colMeans(item_mat_s1)
sigma <- cov(item_mat_s1)
dat_s1$mahal_d2 <- mahalanobis(item_mat_s1, center = mu, cov = sigma)

# Critical chi-sq value (df = 35, alpha = 0.001)
chi2_crit <- qchisq(1 - 0.001, df = 35)
cat(glue("  Critical χ²(35, α=.001) = {round(chi2_crit, 2)}"), "\n")

n_before2 <- nrow(dat_s1)
dat_s2 <- dat_s1 %>% filter(mahal_d2 <= chi2_crit)
n_md_excl <- n_before2 - nrow(dat_s2)
cat(glue("  Mahalanobis D² outliers excluded: n = {n_md_excl} ({round(n_md_excl/n_before2*100,1)}%)"), "\n")
cat(glue("  Remaining: N = {nrow(dat_s2)}"), "\n\n")

# --- Stage 3: Goldbricker Redundancy Check ----------------------------------
cat("──── STAGE 3: Goldbricker Redundancy Check ────\n")

gb_result <- goldbricker(dat_s2[, all_items], p = 0.05, method = "hittner2003",
                         threshold = 0.25, corMin = 0.50, progressbar = FALSE)

# Extract flagged pairs
gb_pairs <- gb_result$proportion_matrix
flagged <- which(gb_pairs > 0.90, arr.ind = TRUE)
flagged <- flagged[flagged[,1] < flagged[,2], , drop = FALSE]  # upper triangle

if (nrow(flagged) == 0) {
  cat("  Goldbricker: Zero item pairs exceeded r > 0.90 threshold.\n")
  cat("  No items excluded for redundancy.\n\n")
} else {
  cat(glue("  Flagged pairs: {nrow(flagged)} — review recommended.\n\n"))
}

# --- Final analytic sample ---------------------------------------------------
dat <- dat_s2
N <- nrow(dat)
cat(strrep("-", 50), "\n")
cat(glue("  FINAL ANALYTIC SAMPLE: N = {N} ({round(N/nrow(raw)*100,1)}% retention)"), "\n")
cat(strrep("-", 50), "\n\n")

# Extract clean item-level data matrix
items_df <- dat[, all_items]

# ──────────────────────────────────────────────────────────────────────────────
# 3.  DESCRIPTIVE STATISTICS
# ──────────────────────────────────────────────────────────────────────────────

cat("──── DESCRIPTIVE STATISTICS ────\n\n")

desc <- describe(items_df)
desc$label <- item_labels[rownames(desc)]
desc$subscale <- subscale_membership[rownames(desc)]
cat("  Item-level summary (selected columns):\n")
print(desc[, c("label", "subscale", "n", "mean", "sd", "skew", "kurtosis")],
      digits = 2)
cat("\n")

# ──────────────────────────────────────────────────────────────────────────────
# 4.  NETWORK ESTIMATION — GGM VIA EBICglasso
# ──────────────────────────────────────────────────────────────────────────────

cat("──── NETWORK ESTIMATION ────\n")

net <- estimateNetwork(items_df,
                       default    = "EBICglasso",
                       corMethod  = "cor_auto",
                       tuning     = 0.5,         # γ = 0.5
                       threshold  = FALSE)

adj_mat  <- getWmat(net)         # Weighted adjacency matrix
n_edges  <- sum(adj_mat[upper.tri(adj_mat)] != 0)
density  <- n_edges / choose(35, 2)
mean_w   <- mean(abs(adj_mat[upper.tri(adj_mat)][adj_mat[upper.tri(adj_mat)] != 0]))
sd_w     <- sd(abs(adj_mat[upper.tri(adj_mat)][adj_mat[upper.tri(adj_mat)] != 0]))

cat(glue("  Edges: {n_edges} / {choose(35,2)} possible (density = {round(density,2)})"), "\n")
cat(glue("  Mean |edge weight|: {round(mean_w,3)} (SD = {round(sd_w,3)})"), "\n\n")

# Visualisation
groups_vis <- list(
  Somatization = which(all_items %in% som_items),
  Depression   = which(all_items %in% dep_items),
  Anxiety      = which(all_items %in% anx_items)
)

pdf("Fig1_network_plot.pdf", width = 12, height = 10)
qgraph(adj_mat,
       layout    = "spring",
       groups    = groups_vis,
       labels    = gsub("Q", "", all_items),
       nodeNames = item_labels[all_items],
       legend.cex = 0.5,
       title     = "GGM (EBICglasso, γ = 0.5)",
       color     = c("#E74C3C", "#3498DB", "#2ECC71"),
       border.width = 1.5,
       vsize     = 6)
dev.off()
cat("  → Saved: Fig1_network_plot.pdf\n\n")

# ──────────────────────────────────────────────────────────────────────────────
# 5.  HYPOTHESIS 1 — STRUCTURAL FUSION (Community Detection + ARI)
# ──────────────────────────────────────────────────────────────────────────────

cat(strrep("=", 78), "\n")
cat("  HYPOTHESIS 1: Structural Fusion\n")
cat(strrep("=", 78), "\n\n")

# Convert adjacency to igraph (absolute weights for community detection)
g <- graph_from_adjacency_matrix(abs(adj_mat), mode = "undirected",
                                 weighted = TRUE, diag = FALSE)

# Walktrap (steps = 4)
wt <- cluster_walktrap(g, steps = 4)
wt_membership <- membership(wt)

# Louvain
lv <- cluster_louvain(g, weights = E(g)$weight)
lv_membership <- membership(lv)

# Modularity
Q_wt <- modularity(wt)
Q_lv <- modularity(lv)

cat(glue("  Walktrap:  {max(wt_membership)} communities, Q = {round(Q_wt, 3)}"), "\n")
cat(glue("  Louvain:   {max(lv_membership)} communities, Q = {round(Q_lv, 3)}"), "\n\n")

# Community assignments table
comm_table <- data.frame(
  Item      = all_items,
  Label     = item_labels[all_items],
  SCL_Sub   = subscale_membership[all_items],
  Walktrap  = wt_membership,
  Louvain   = lv_membership,
  stringsAsFactors = FALSE
)
cat("  Community assignments:\n")
print(comm_table, row.names = FALSE)
cat("\n")

# Agreement between Walktrap and Louvain
ari_wt_lv <- adjustedRandIndex(wt_membership, lv_membership)
cat(glue("  ARI (Walktrap vs Louvain): {round(ari_wt_lv, 3)}"), "\n")

# ARI: Empirical communities vs SCL-90-R 3-subscale structure
scl_membership_numeric <- as.numeric(factor(subscale_membership[all_items],
                                            levels = c("Somatization", "Depression", "Anxiety")))
ari_wt_scl <- adjustedRandIndex(wt_membership, scl_membership_numeric)
ari_lv_scl <- adjustedRandIndex(lv_membership, scl_membership_numeric)

cat(glue("  ARI (Walktrap vs SCL-90-R subscales): {round(ari_wt_scl, 3)}"), "\n")
cat(glue("  ARI (Louvain  vs SCL-90-R subscales): {round(ari_lv_scl, 3)}"), "\n\n")

# ARI: Empirical communities vs hypothesised 2-domain Somatic-Affective model
hyp_2dom <- as.numeric(factor(bridge_communities[all_items]))
ari_wt_2dom <- adjustedRandIndex(wt_membership, hyp_2dom)
ari_lv_2dom <- adjustedRandIndex(lv_membership, hyp_2dom)

cat(glue("  ARI (Walktrap vs 2-domain Somatic-Affective): {round(ari_wt_2dom, 3)}"), "\n")
cat(glue("  ARI (Louvain  vs 2-domain Somatic-Affective): {round(ari_lv_2dom, 3)}"), "\n\n")

# H1 decision
h1_support <- (ari_wt_scl < 0.50) & (ari_lv_scl < 0.50)
cat(glue("  ▶ H1 DECISION: ARI with SCL subscales < 0.50? {ifelse(h1_support, 'YES — Structural Fusion SUPPORTED', 'NO — Falsified')}"), "\n\n")

# ──────────────────────────────────────────────────────────────────────────────
# 6.  HYPOTHESIS 2 — BRIDGE ASYMMETRY (Bridge Expected Influence)
# ──────────────────────────────────────────────────────────────────────────────

cat(strrep("=", 78), "\n")
cat("  HYPOTHESIS 2: Bridge Asymmetry\n")
cat(strrep("=", 78), "\n\n")

# Define a priori communities for bridge analysis (Somatic vs Affective)
bridge_comm_vec <- ifelse(all_items %in% som_items, 1, 2)
names(bridge_comm_vec) <- all_items

# Calculate Bridge Expected Influence (1-step)
bei <- bridge(adj_mat, communities = bridge_comm_vec, useCommunities = "all")
bei_values <- bei$`Bridge Expected Influence (1-step)`

# Attach domain labels
bei_df <- data.frame(
  Item     = all_items,
  Label    = item_labels[all_items],
  Domain   = bridge_communities[all_items],
  BEI      = bei_values,
  stringsAsFactors = FALSE
) %>% arrange(desc(abs(BEI)))

cat("  Bridge Expected Influence (sorted by |BEI|):\n")
print(bei_df, row.names = FALSE, digits = 3)
cat("\n")

# Group comparison: Somatic vs Affective BEI
som_bei <- bei_df$BEI[bei_df$Domain == "Somatic"]
aff_bei <- bei_df$BEI[bei_df$Domain == "Affective"]

t_result <- t.test(som_bei, aff_bei, var.equal = FALSE)
d_result <- cohens_d(som_bei, aff_bei)

cat(glue("  Mean BEI — Somatic: {round(mean(som_bei), 4)}, Affective: {round(mean(aff_bei), 4)}"), "\n")
cat(glue("  t({round(t_result$parameter, 1)}) = {round(t_result$statistic, 3)}, p = {format.pval(t_result$p.value, digits = 3)}"), "\n")
cat(glue("  Cohen's d = {round(d_result$Cohens_d, 3)} [{round(d_result$CI_low, 3)}, {round(d_result$CI_high, 3)}]"), "\n\n")

h2_support <- (d_result$Cohens_d >= 0.80) & (t_result$p.value < 0.01)
cat(glue("  ▶ H2 DECISION: d ≥ 0.80 & p < .01? {ifelse(h2_support, 'YES — Bridge Asymmetry SUPPORTED', 'NO — Falsified')}"), "\n\n")

# Bootstrap BEI CIs (5,000 iterations)
cat("  Running bootstrapped BEI CIs (5000 iterations)... ")
n_boot <- 5000
boot_bei_mat <- matrix(NA, nrow = n_boot, ncol = 35)
colnames(boot_bei_mat) <- all_items

for (b in seq_len(n_boot)) {
  idx <- sample(nrow(items_df), replace = TRUE)
  boot_net <- suppressWarnings(
    estimateNetwork(items_df[idx, ], default = "EBICglasso",
                    corMethod = "cor_auto", tuning = 0.5, threshold = FALSE,
                    verbose = FALSE)
  )
  boot_adj <- getWmat(boot_net)
  boot_bridge <- bridge(boot_adj, communities = bridge_comm_vec, useCommunities = "all")
  boot_bei_mat[b, ] <- boot_bridge$`Bridge Expected Influence (1-step)`
  if (b %% 1000 == 0) cat(glue("{b}.."))
}
cat(" done.\n")

bei_ci <- data.frame(
  Item  = all_items,
  Label = item_labels[all_items],
  BEI   = bei_values,
  CI_lo = apply(boot_bei_mat, 2, quantile, 0.025),
  CI_hi = apply(boot_bei_mat, 2, quantile, 0.975)
)
cat("  Bootstrapped 95% CIs for BEI:\n")
print(bei_ci %>% arrange(desc(abs(BEI))), row.names = FALSE, digits = 3)
cat("\n")

# Bridge plot
pdf("Fig2_bridge_bei.pdf", width = 10, height = 7)
par(mar = c(5, 10, 4, 2))
ord <- order(bei_values)
cols <- ifelse(all_items[ord] %in% som_items, "#E74C3C", "#3498DB")
barplot(bei_values[ord], names.arg = item_labels[all_items[ord]],
        horiz = TRUE, las = 1, col = cols, border = NA,
        main = "Bridge Expected Influence (1-step)",
        xlab = "BEI", cex.names = 0.6)
legend("bottomright", legend = c("Somatic", "Affective"),
       fill = c("#E74C3C", "#3498DB"), bty = "n")
dev.off()
cat("  → Saved: Fig2_bridge_bei.pdf\n\n")

# ──────────────────────────────────────────────────────────────────────────────
# 7.  HYPOTHESIS 3 — MECHANISTIC ESSENTIALITY (In Silico Node Knockout)
# ──────────────────────────────────────────────────────────────────────────────

cat(strrep("=", 78), "\n")
cat("  HYPOTHESIS 3: Mechanistic Essentiality (Node Knockout Simulation)\n")
cat(strrep("=", 78), "\n\n")

# --- Custom knockout function ------------------------------------------------
node_knockout <- function(data, items, target_node, gamma = 0.5) {
  #' Removes target_node from data, re-estimates GGM, returns ΔConnectivity.
  #'

  #' @param data       Data frame of item responses.
  #' @param items      Character vector of all item names.
  #' @param target_node Character: name of node to knock out.
  #' @param gamma      Numeric: EBIC tuning parameter.
  #'
  #' @return List with baseline_S, knockout_S, delta_connectivity (%).

  # Baseline connectivity
  net_base <- estimateNetwork(data[, items], default = "EBICglasso",
                              corMethod = "cor_auto", tuning = gamma,
                              threshold = FALSE, verbose = FALSE)
  adj_base <- getWmat(net_base)
  S_base <- sum(abs(adj_base[upper.tri(adj_base)]))

  # Remove target node
  remaining <- setdiff(items, target_node)
  net_ko <- estimateNetwork(data[, remaining], default = "EBICglasso",
                            corMethod = "cor_auto", tuning = gamma,
                            threshold = FALSE, verbose = FALSE)
  adj_ko <- getWmat(net_ko)
  S_ko <- sum(abs(adj_ko[upper.tri(adj_ko)]))

  delta <- (S_base - S_ko) / S_base * 100

  list(
    target_node      = target_node,
    baseline_S       = S_base,
    knockout_S       = S_ko,
    delta_connectivity = delta
  )
}

# --- Identify top-6 somatic and centrality-matched affective nodes -----------
cat("  Identifying top-6 BEI somatic nodes and centrality-matched affective nodes...\n")

# Top 6 somatic by BEI
top6_som <- bei_df %>%
  filter(Domain == "Somatic") %>%
  arrange(desc(abs(BEI))) %>%
  slice_head(n = 6) %>%
  pull(Item)

cat(glue("  Top-6 Somatic (BEI): {paste(top6_som, collapse = ', ')}"), "\n")

# Expected Influence for matching
ei_vals <- centrality_auto(adj_mat)$node.centrality$ExpectedInfluence
names(ei_vals) <- all_items

# For each somatic node, find closest-EI affective node (±0.10)
som_ei <- ei_vals[top6_som]
aff_candidates <- bei_df %>% filter(Domain == "Affective") %>% pull(Item)

matched_aff <- character(6)
used <- c()
for (i in seq_along(top6_som)) {
  diffs <- abs(ei_vals[aff_candidates] - som_ei[i])
  diffs[aff_candidates %in% used] <- Inf
  best <- names(which.min(diffs))
  matched_aff[i] <- best
  used <- c(used, best)
  cat(glue("    {top6_som[i]} (EI={round(som_ei[i],3)}) ↔ {best} (EI={round(ei_vals[best],3)}, Δ={round(diffs[best],3)})"), "\n")
}
cat("\n")

# --- Run knockout loop -------------------------------------------------------
cat("  Running knockout simulations (12 nodes)...\n")

knockout_targets <- c(top6_som, matched_aff)
knockout_labels  <- c(rep("Somatic", 6), rep("Affective", 6))

ko_results <- vector("list", length(knockout_targets))
for (i in seq_along(knockout_targets)) {
  node <- knockout_targets[i]
  cat(glue("    [{i}/12] Knocking out {node} ({item_labels[node]})..."))
  ko_results[[i]] <- node_knockout(items_df, all_items, node, gamma = 0.5)
  cat(glue(" ΔConnectivity = {round(ko_results[[i]]$delta_connectivity, 2)}%"), "\n")
}
cat("\n")

# Compile results
ko_df <- do.call(rbind, lapply(seq_along(ko_results), function(i) {
  data.frame(
    Node   = ko_results[[i]]$target_node,
    Label  = item_labels[ko_results[[i]]$target_node],
    Domain = knockout_labels[i],
    S_base = ko_results[[i]]$baseline_S,
    S_ko   = ko_results[[i]]$knockout_S,
    Delta  = ko_results[[i]]$delta_connectivity,
    stringsAsFactors = FALSE
  )
}))

cat("  Knockout results:\n")
print(ko_df, row.names = FALSE, digits = 3)
cat("\n")

# Statistical comparison: Somatic vs Affective ΔConnectivity
som_delta <- ko_df$Delta[ko_df$Domain == "Somatic"]
aff_delta <- ko_df$Delta[ko_df$Domain == "Affective"]

# Paired permutation test (5,000 iterations)
obs_diff <- mean(som_delta) - mean(aff_delta)
n_perm <- 5000
perm_diffs <- numeric(n_perm)
all_deltas <- c(som_delta, aff_delta)
for (p in seq_len(n_perm)) {
  perm_idx <- sample(12, 6)
  perm_diffs[p] <- mean(all_deltas[perm_idx]) - mean(all_deltas[-perm_idx])
}
perm_p <- mean(abs(perm_diffs) >= abs(obs_diff))

d_ko <- cohens_d(som_delta, aff_delta)

cat(glue("  Mean ΔConnectivity — Somatic: {round(mean(som_delta), 2)}%, Affective: {round(mean(aff_delta), 2)}%"), "\n")
cat(glue("  Difference: {round(obs_diff, 2)} percentage points"), "\n")
cat(glue("  Permutation p = {round(perm_p, 4)}"), "\n")
cat(glue("  Cohen's d = {round(d_ko$Cohens_d, 3)}"), "\n\n")

# Null distribution from 1,000 random 6-node samples
cat("  Computing null distribution (1,000 random 6-node knockouts)...\n")
null_deltas <- numeric(1000)
for (r in seq_len(1000)) {
  rand_nodes <- sample(all_items, 6)
  rand_ko <- sapply(rand_nodes, function(nd) {
    res <- node_knockout(items_df, all_items, nd, gamma = 0.5)
    res$delta_connectivity
  })
  null_deltas[r] <- mean(rand_ko)
  if (r %% 100 == 0) cat(glue("    {r}/1000.."), "\n")
}

null_p_som <- mean(null_deltas >= mean(som_delta))
null_p_aff <- mean(null_deltas >= mean(aff_delta))
cat(glue("\n  Null mean ΔConnectivity: {round(mean(null_deltas), 2)}% (SD = {round(sd(null_deltas), 2)}%)"), "\n")
cat(glue("  Somatic vs null p = {round(null_p_som, 4)}"), "\n")
cat(glue("  Affective vs null p = {round(null_p_aff, 4)}"), "\n\n")

h3_support <- (obs_diff >= 10) & (perm_p < 0.01) & (d_ko$Cohens_d >= 1.00)
cat(glue("  ▶ H3 DECISION: Δ ≥ 10 pp, p < .01, d ≥ 1.0? {ifelse(h3_support, 'YES — Mechanistic Essentiality SUPPORTED', 'NO — Falsified')}"), "\n\n")

# Knockout visualisation
pdf("Fig3_knockout_barplot.pdf", width = 9, height = 6)
ggplot(ko_df, aes(x = reorder(Label, Delta), y = Delta, fill = Domain)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("Somatic" = "#E74C3C", "Affective" = "#3498DB")) +
  labs(title = "In Silico Node Knockout: ΔConnectivity (%)",
       x = NULL, y = "ΔConnectivity (%)") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")
dev.off()
cat("  → Saved: Fig3_knockout_barplot.pdf\n\n")

# ──────────────────────────────────────────────────────────────────────────────
# 8.  HYPOTHESIS 4 — SENTINEL SHORT-FORM VALIDATION
# ──────────────────────────────────────────────────────────────────────────────

cat(strrep("=", 78), "\n")
cat("  HYPOTHESIS 4: Sentinel Short-Form (Information-Theoretic Compression)\n")
cat(strrep("=", 78), "\n\n")

# Step 1: Identify 6 highest-BEI somatic nodes (already computed)
sentinel_items <- bei_df %>%
  filter(Domain == "Somatic") %>%
  arrange(desc(abs(BEI))) %>%
  slice_head(n = 6) %>%
  pull(Item)

cat(glue("  Sentinel items: {paste(sentinel_items, collapse = ', ')}"), "\n")
cat(glue("  Labels: {paste(item_labels[sentinel_items], collapse = ', ')}"), "\n\n")

# Step 2: Derivation / Validation split (60/40)
set.seed(2024)
deriv_idx <- sample(nrow(items_df), size = round(0.60 * nrow(items_df)))
valid_idx <- setdiff(seq_len(nrow(items_df)), deriv_idx)

deriv_data <- items_df[deriv_idx, ]
valid_data <- items_df[valid_idx, ]

cat(glue("  Derivation: n = {nrow(deriv_data)}, Validation: n = {nrow(valid_data)}"), "\n\n")

# Full-scale internalizing score (sum of 35 items)
deriv_data$full_score <- rowSums(deriv_data[, all_items])
valid_data$full_score <- rowSums(valid_data[, all_items])

# Sentinel sum score
deriv_data$sentinel_score <- rowSums(deriv_data[, sentinel_items])
valid_data$sentinel_score <- rowSums(valid_data[, sentinel_items])

# Step 3: Criterion validation — R²
lm_sentinel <- lm(full_score ~ sentinel_score, data = valid_data)
r2_sentinel <- summary(lm_sentinel)$r.squared

cat(glue("  Sentinel R² (validation): {round(r2_sentinel, 4)}"), "\n")

# ROC analysis — clinical discrimination at T-score ≥ 63
# T-score conversion: T = 50 + 10 * (X - M) / SD
full_mean <- mean(c(deriv_data$full_score, valid_data$full_score))
full_sd   <- sd(c(deriv_data$full_score, valid_data$full_score))
valid_data$t_score <- 50 + 10 * (valid_data$full_score - full_mean) / full_sd
valid_data$clinical_case <- ifelse(valid_data$t_score >= 63, 1, 0)

roc_sentinel <- roc(valid_data$clinical_case, valid_data$sentinel_score,
                    quiet = TRUE)
auc_sentinel <- auc(roc_sentinel)
coords_sentinel <- coords(roc_sentinel, "best", ret = c("threshold", "sensitivity",
                                                          "specificity", "ppv", "npv"))

cat(glue("  AUC = {round(auc_sentinel, 3)}"), "\n")
cat(glue("  Optimal threshold = {round(coords_sentinel$threshold, 2)} "),
    glue("(Sens = {round(coords_sentinel$sensitivity, 3)}, "),
    glue("Spec = {round(coords_sentinel$specificity, 3)})"), "\n\n")

# Step 4: Comparison with four alternative short-forms
cat("  Comparing with alternative short-forms...\n")

# (a) Random 6-item somatic samples (1,000 permutations)
set.seed(2024)
random_r2 <- replicate(1000, {
  rand_items <- sample(som_items, 6)
  valid_data$rand_score <- rowSums(valid_data[, rand_items])
  summary(lm(full_score ~ rand_score, data = valid_data))$r.squared
})

# (b) 6 highest centrality (non-bridge, i.e. highest EI) somatic nodes
top6_ei_som <- names(sort(ei_vals[som_items], decreasing = TRUE))[1:6]
valid_data$ei_score <- rowSums(valid_data[, top6_ei_som])
r2_ei <- summary(lm(full_score ~ ei_score, data = valid_data))$r.squared

# (c) Balanced: 2 somatic + 2 depression + 2 anxiety (highest BEI each)
top2_som <- bei_df %>% filter(Domain == "Somatic") %>%
  arrange(desc(abs(BEI))) %>% slice_head(n = 2) %>% pull(Item)
# For affective, split into Dep and Anx
bei_dep <- bei_df %>% filter(Item %in% dep_items) %>%
  arrange(desc(abs(BEI))) %>% slice_head(n = 2) %>% pull(Item)
bei_anx <- bei_df %>% filter(Item %in% anx_items) %>%
  arrange(desc(abs(BEI))) %>% slice_head(n = 2) %>% pull(Item)
balanced_items <- c(top2_som, bei_dep, bei_anx)
valid_data$balanced_score <- rowSums(valid_data[, balanced_items])
r2_balanced <- summary(lm(full_score ~ balanced_score, data = valid_data))$r.squared

# (d) 6 highest-EI affective symptoms
top6_ei_aff <- names(sort(ei_vals[c(dep_items, anx_items)], decreasing = TRUE))[1:6]
valid_data$aff_score <- rowSums(valid_data[, top6_ei_aff])
r2_aff <- summary(lm(full_score ~ aff_score, data = valid_data))$r.squared

# Summary table
comparison_table <- data.frame(
  ShortForm     = c("Sentinel (Bridge Somatic)", "Random Somatic (mean)",
                    "Top Centrality Somatic", "Balanced (2+2+2)",
                    "Top Affective"),
  R2            = c(r2_sentinel, mean(random_r2), r2_ei, r2_balanced, r2_aff),
  Delta_R2_vs_Sentinel = c(NA, r2_sentinel - mean(random_r2),
                            r2_sentinel - r2_ei, r2_sentinel - r2_balanced,
                            r2_sentinel - r2_aff)
)

cat("\n  Short-form comparison (Validation sample):\n")
print(comparison_table, row.names = FALSE, digits = 4)
cat("\n")

# Bootstrapped ΔR² significance (5,000 iterations, Bonferroni α = .01/4 = .0025)
cat("  Running bootstrapped ΔR² tests (5,000 iterations)...\n")
n_boot_r2 <- 5000
boot_delta_ei  <- numeric(n_boot_r2)
boot_delta_bal <- numeric(n_boot_r2)
boot_delta_aff <- numeric(n_boot_r2)

for (b in seq_len(n_boot_r2)) {
  idx <- sample(nrow(valid_data), replace = TRUE)
  bsample <- valid_data[idx, ]

  r2_sent_b <- summary(lm(full_score ~ sentinel_score, data = bsample))$r.squared
  r2_ei_b   <- summary(lm(full_score ~ ei_score, data = bsample))$r.squared
  r2_bal_b  <- summary(lm(full_score ~ balanced_score, data = bsample))$r.squared
  r2_aff_b  <- summary(lm(full_score ~ aff_score, data = bsample))$r.squared

  boot_delta_ei[b]  <- r2_sent_b - r2_ei_b
  boot_delta_bal[b] <- r2_sent_b - r2_bal_b
  boot_delta_aff[b] <- r2_sent_b - r2_aff_b
}

p_vs_ei  <- mean(boot_delta_ei  <= 0) * 4  # Bonferroni
p_vs_bal <- mean(boot_delta_bal <= 0) * 4
p_vs_aff <- mean(boot_delta_aff <= 0) * 4

cat(glue("  ΔR² Sentinel vs Centrality: {round(r2_sentinel - r2_ei, 4)}, p = {round(min(p_vs_ei,1), 4)}"), "\n")
cat(glue("  ΔR² Sentinel vs Balanced:   {round(r2_sentinel - r2_balanced, 4)}, p = {round(min(p_vs_bal,1), 4)}"), "\n")
cat(glue("  ΔR² Sentinel vs Affective:  {round(r2_sentinel - r2_aff, 4)}, p = {round(min(p_vs_aff,1), 4)}"), "\n\n")

h4_support <- (r2_sentinel >= 0.80) &
  (r2_sentinel - r2_ei >= 0.10) &
  (r2_sentinel - r2_balanced >= 0.10) &
  (r2_sentinel - r2_aff >= 0.10)
cat(glue("  ▶ H4 DECISION: R² ≥ 0.80 & ΔR² ≥ 0.10? {ifelse(h4_support, 'YES — Sentinel SUPPORTED', 'NO — Falsified')}"), "\n\n")

# ROC plot
pdf("Fig4_sentinel_roc.pdf", width = 7, height = 7)
plot(roc_sentinel, main = "Sentinel Short-Form: ROC Curve",
     col = "#E74C3C", lwd = 2, print.auc = TRUE, print.auc.x = 0.4,
     legacy.axes = TRUE)
dev.off()
cat("  → Saved: Fig4_sentinel_roc.pdf\n\n")

# ──────────────────────────────────────────────────────────────────────────────
# 9.  HYPOTHESIS 5 — TEMPORAL INVARIANCE (NCT)
# ──────────────────────────────────────────────────────────────────────────────

cat(strrep("=", 78), "\n")
cat("  HYPOTHESIS 5: Temporal Invariance (Network Comparison Test)\n")
cat(strrep("=", 78), "\n\n")

if (all(is.na(dat$cohort_year))) {
  cat("  ⚠ Skipping H5: No cohort_year variable available.\n\n")
} else {

  # Define three pandemic-aligned periods
  # Adapt these year values to match your data coding
  period_map <- list(
    "Pre-pandemic"  = c(2019, 2020),
    "Acute"         = c(2021, 2022),
    "Post-acute"    = c(2023, 2024)
  )

  dat$period <- case_when(
    dat$cohort_year %in% period_map$`Pre-pandemic` ~ "Pre-pandemic",
    dat$cohort_year %in% period_map$Acute          ~ "Acute",
    dat$cohort_year %in% period_map$`Post-acute`   ~ "Post-acute",
    TRUE ~ NA_character_
  )

  period_ns <- table(dat$period)
  cat("  Period sample sizes:\n")
  print(period_ns)
  cat("\n")

  # Prepare data for each period
  dat_pre  <- dat %>% filter(period == "Pre-pandemic") %>% select(all_of(all_items))
  dat_acu  <- dat %>% filter(period == "Acute")        %>% select(all_of(all_items))
  dat_post <- dat %>% filter(period == "Post-acute")   %>% select(all_of(all_items))

  # --- NCT: Pre-pandemic vs Acute ---
  cat("  Running NCT: Pre-pandemic vs Acute (2500 permutations)...\n")
  nct_pre_acu <- NCT(dat_pre, dat_acu,
                     it = 2500,
                     binary.data = FALSE,
                     paired = FALSE,
                     weighted = TRUE,
                     test.edges = TRUE,
                     edges = "all",
                     progressbar = TRUE,
                     p.adjust.methods = "BH",
                     gamma = 0.5)

  cat(glue("  Network structure: M = {round(nct_pre_acu$nw.invariance.real, 4)}, p = {round(nct_pre_acu$nw.invariance.pval, 4)}"), "\n")
  cat(glue("  Global strength:   S = {round(nct_pre_acu$glstrinv.real, 4)}, p = {round(nct_pre_acu$glstrinv.pval, 4)}"), "\n")
  n_sig_edges_1 <- sum(nct_pre_acu$einv.pvals$`p-value` < 0.05, na.rm = TRUE)
  cat(glue("  Significant edge differences (FDR < .05): {n_sig_edges_1}"), "\n\n")

  # --- NCT: Pre-pandemic vs Post-acute ---
  cat("  Running NCT: Pre-pandemic vs Post-acute (2500 permutations)...\n")
  nct_pre_post <- NCT(dat_pre, dat_post,
                      it = 2500,
                      binary.data = FALSE,
                      paired = FALSE,
                      weighted = TRUE,
                      test.edges = TRUE,
                      edges = "all",
                      progressbar = TRUE,
                      p.adjust.methods = "BH",
                      gamma = 0.5)

  cat(glue("  Network structure: M = {round(nct_pre_post$nw.invariance.real, 4)}, p = {round(nct_pre_post$nw.invariance.pval, 4)}"), "\n")
  cat(glue("  Global strength:   S = {round(nct_pre_post$glstrinv.real, 4)}, p = {round(nct_pre_post$glstrinv.pval, 4)}"), "\n")
  n_sig_edges_2 <- sum(nct_pre_post$einv.pvals$`p-value` < 0.05, na.rm = TRUE)
  cat(glue("  Significant edge differences (FDR < .05): {n_sig_edges_2}"), "\n\n")

  # --- NCT: Acute vs Post-acute ---
  cat("  Running NCT: Acute vs Post-acute (2500 permutations)...\n")
  nct_acu_post <- NCT(dat_acu, dat_post,
                      it = 2500,
                      binary.data = FALSE,
                      paired = FALSE,
                      weighted = TRUE,
                      test.edges = TRUE,
                      edges = "all",
                      progressbar = TRUE,
                      p.adjust.methods = "BH",
                      gamma = 0.5)

  cat(glue("  Network structure: M = {round(nct_acu_post$nw.invariance.real, 4)}, p = {round(nct_acu_post$nw.invariance.pval, 4)}"), "\n")
  cat(glue("  Global strength:   S = {round(nct_acu_post$glstrinv.real, 4)}, p = {round(nct_acu_post$glstrinv.pval, 4)}"), "\n")
  n_sig_edges_3 <- sum(nct_acu_post$einv.pvals$`p-value` < 0.05, na.rm = TRUE)
  cat(glue("  Significant edge differences (FDR < .05): {n_sig_edges_3}"), "\n\n")

  # --- Cross-period community stability ---
  cat("  Cross-period community detection (Walktrap)...\n")

  net_pre  <- estimateNetwork(dat_pre,  default = "EBICglasso", corMethod = "cor_auto",
                              tuning = 0.5, verbose = FALSE)
  net_acu  <- estimateNetwork(dat_acu,  default = "EBICglasso", corMethod = "cor_auto",
                              tuning = 0.5, verbose = FALSE)
  net_post <- estimateNetwork(dat_post, default = "EBICglasso", corMethod = "cor_auto",
                              tuning = 0.5, verbose = FALSE)

  g_pre  <- graph_from_adjacency_matrix(abs(getWmat(net_pre)),  mode = "undirected",
                                        weighted = TRUE, diag = FALSE)
  g_acu  <- graph_from_adjacency_matrix(abs(getWmat(net_acu)),  mode = "undirected",
                                        weighted = TRUE, diag = FALSE)
  g_post <- graph_from_adjacency_matrix(abs(getWmat(net_post)), mode = "undirected",
                                        weighted = TRUE, diag = FALSE)

  wt_pre  <- membership(cluster_walktrap(g_pre,  steps = 4))
  wt_acu  <- membership(cluster_walktrap(g_acu,  steps = 4))
  wt_post <- membership(cluster_walktrap(g_post, steps = 4))

  ari_pre_acu  <- adjustedRandIndex(wt_pre, wt_acu)
  ari_pre_post <- adjustedRandIndex(wt_pre, wt_post)
  ari_acu_post <- adjustedRandIndex(wt_acu, wt_post)

  cat(glue("  ARI (Pre vs Acute):     {round(ari_pre_acu, 3)}"), "\n")
  cat(glue("  ARI (Pre vs Post):      {round(ari_pre_post, 3)}"), "\n")
  cat(glue("  ARI (Acute vs Post):    {round(ari_acu_post, 3)}"), "\n\n")

  # H5 decision
  struct_invariant <- (nct_pre_acu$nw.invariance.pval > 0.05) &
    (nct_pre_post$nw.invariance.pval > 0.05) &
    (nct_acu_post$nw.invariance.pval > 0.05)
  comm_stable <- (ari_pre_acu > 0.80) & (ari_pre_post > 0.80) & (ari_acu_post > 0.80)

  h5_support <- struct_invariant & comm_stable
  cat(glue("  ▶ H5 DECISION: Network invariance (all p > .05) & community ARI > 0.80?"), "\n")
  cat(glue("    Structure invariant: {struct_invariant}"), "\n")
  cat(glue("    Community stable:    {comm_stable}"), "\n")
  cat(glue("    {ifelse(h5_support, 'YES — Temporal Invariance SUPPORTED', 'NO — Falsified')}"), "\n\n")
}

# ──────────────────────────────────────────────────────────────────────────────
# 10.  COMPREHENSIVE SUMMARY REPORT
# ──────────────────────────────────────────────────────────────────────────────

cat("\n", strrep("█", 78), "\n")
cat("  ANALYSIS COMPLETE — SUMMARY REPORT\n")
cat(strrep("█", 78), "\n\n")

cat(glue("  Final N = {N} ({round(N/nrow(raw)*100,1)}% retention from N = {nrow(raw)})"), "\n")
cat(glue("  Network: {n_edges} edges, density = {round(density, 2)}, mean |w| = {round(mean_w, 3)}"), "\n\n")

cat("  ┌──────────────────────────────────────────────────────────────────────┐\n")
cat("  │  HYPOTHESIS                          │  RESULT                      │\n")
cat("  ├──────────────────────────────────────────────────────────────────────┤\n")
cat(glue("  │  H1: Structural Fusion               │  ARI(WT)={round(ari_wt_scl,2)}, ARI(LV)={round(ari_lv_scl,2)}"), "\n")
cat(glue("  │      (ARI < 0.50 → supported)        │  {ifelse(h1_support, 'SUPPORTED', 'FALSIFIED')}"), "\n")
cat(glue("  │  H2: Bridge Asymmetry                │  d = {round(d_result$Cohens_d,2)}, p = {format.pval(t_result$p.value, digits=3)}"), "\n")
cat(glue("  │      (d ≥ 0.80, p < .01)             │  {ifelse(h2_support, 'SUPPORTED', 'FALSIFIED')}"), "\n")
cat(glue("  │  H3: Mechanistic Essentiality        │  Δ = {round(obs_diff,2)}pp, d = {round(d_ko$Cohens_d,2)}"), "\n")
cat(glue("  │      (Δ ≥ 10pp, d ≥ 1.00, p < .01)  │  {ifelse(h3_support, 'SUPPORTED', 'FALSIFIED')}"), "\n")
cat(glue("  │  H4: Sentinel Short-Form             │  R² = {round(r2_sentinel,3)}"), "\n")
cat(glue("  │      (R² ≥ 0.80, ΔR² ≥ 0.10)        │  {ifelse(h4_support, 'SUPPORTED', 'FALSIFIED')}"), "\n")

if (!all(is.na(dat$cohort_year))) {
  cat(glue("  │  H5: Temporal Invariance             │  {ifelse(h5_support, 'SUPPORTED', 'FALSIFIED')}"), "\n")
} else {
  cat("  │  H5: Temporal Invariance             │  SKIPPED (no year data)      │\n")
}

cat("  └──────────────────────────────────────────────────────────────────────┘\n\n")

cat("  Output files:\n")
cat("    • Fig1_network_plot.pdf\n")
cat("    • Fig2_bridge_bei.pdf\n")
cat("    • Fig3_knockout_barplot.pdf\n")
cat("    • Fig4_sentinel_roc.pdf\n\n")

cat("  Session info:\n")
cat(glue("    R version: {R.version.string}"), "\n")
cat(glue("    Date:      {Sys.Date()}"), "\n")
cat(glue("    Platform:  {R.version$platform}"), "\n\n")

cat(strrep("=", 78), "\n")
cat("  END OF ANALYSIS\n")
cat(strrep("=", 78), "\n")

