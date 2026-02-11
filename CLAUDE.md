# CLAUDE.md

## Project Overview

This repository contains a **Network Analysis and In-Silico Intervention Simulation** research project studying psychological symptom networks in a Chinese university population (N = 47,860). The study uses network psychometrics to analyze interconnections between somatization, depression, and anxiety symptoms from the SCL-90-R, testing a "Somatic-Affective fusion" hypothesis in internalizing psychopathology.

**Core research questions** are organized around five hypotheses:
- **H1 (Structural Fusion):** Empirical symptom communities differ from traditional SCL-90-R subscales
- **H2 (Bridge Asymmetry):** Somatic symptoms have higher bridge centrality than affective symptoms
- **H3 (Mechanistic Essentiality):** High-bridge somatic nodes cause disproportionate network destabilization when removed
- **H4 (Information-Theoretic Compression):** A 6-item bridge-based screening tool captures core internalizing distress
- **H5 (Temporal Invariance):** Network structure remains stable across pandemic-aligned time periods (2019-2025)

## Repository Structure

```
A-Network-Analysis-and-In-Silico-Intervention-Simulation/
├── CLAUDE.md          # This file - project guide for AI assistants
└── NCT                # Methods documentation (research manuscript section)
```

The repository is in early-stage development. The `NCT` file contains the full methods section describing the analytical pipeline. Analysis scripts and data are referenced but not yet committed (placeholders: `[REPOSITORY]`, `[OSF LINK]`).

## Technology Stack

- **Language:** R 4.3.1
- **Key R packages:**
  - `bootnet` 1.5 — Network estimation (GGM via EBICglasso)
  - `qgraph` 1.9.5 — Network visualization and analysis
  - `igraph` 1.5.0 — Community detection (Walktrap, Louvain algorithms)
  - `NetworkComparisonTest` 2.2.2 — Temporal invariance permutation testing
  - `networktools` 1.5.0 — Bridge centrality measures

## Analytical Pipeline

The analysis follows a sequential pipeline as documented in the `NCT` file:

1. **Data Quality Assurance** — Three-stage filtering (longstring analysis, Mahalanobis distance outliers, Goldbricker redundancy)
2. **Network Estimation** — Gaussian Graphical Models via EBICglasso (gamma = 0.5)
3. **Community Detection** — Walktrap and Louvain algorithms with ARI comparison
4. **Bridge Centrality Analysis** — Bridge expected influence (BEI) with bootstrapped CIs (5,000 iterations)
5. **In-Silico Node Knockout** — Computational removal experiments measuring network destabilization (delta connectivity)
6. **Sentinel Short-Form Development** — 6-item bridge-based screening tool (Q56, Q53, Q52, Q58, Q54, Q51)
7. **Temporal Invariance Testing** — NCT across pre-pandemic, acute pandemic, and post-acute periods

## Key Parameters and Thresholds

| Parameter | Value | Context |
|-----------|-------|---------|
| EBICglasso gamma | 0.5 (default; sensitivity: 0.25, 0.75) | Network sparsity tuning |
| Network density | 0.58 (342 non-zero edges of 595 possible) | Estimated network |
| Community strength | Q > 0.30 | Modularity threshold |
| ARI correspondence | < 0.50 weak, > 0.70 strong, > 0.80 stable | Community comparison |
| Bootstrap iterations | 5,000 (BEI CIs), 1,000 (null distributions) | Resampling |
| NCT permutations | 2,500 | Temporal invariance testing |
| Longstring exclusion | >= 15 consecutive identical responses | Data quality |
| Mahalanobis exclusion | chi-squared critical (df=35, alpha=.001) | Outlier detection |
| Clinical threshold | T-score >= 63 | ROC discrimination |
| Short-form success | R-squared >= 0.80 and delta-R-squared >= 0.10 (p < .01) | Sentinel validation |
| Alpha level | 0.05 (FDR-corrected where applicable) | Statistical significance |

## Data Description

- **Source:** Sichuan Preschool Education College annual mental health screening (2019-2025)
- **Instrument:** Chinese-validated SCL-90-R (Wang, 1984), 35 items from 3 subscales
- **Sample:** N = 47,860 (62.3% female; mean age = 19.7, SD = 1.4), 95.7% retention after quality filtering
- **Subscales:** Somatization (12 items), Depression (13 items), Anxiety (10 items)
- **Temporal cohorts for invariance testing:**
  - Pre-pandemic: 2019-2020 (n = 17,125)
  - Acute pandemic: 2021-2022 (n = 15,234)
  - Post-acute: 2023-2024 (n = 15,501)
- **Data is de-identified** per IRB approval and Declaration of Helsinki

## Development Conventions

### Code Style (R)
- Use R 4.3.1 for all analyses
- Follow tidyverse style guide for R code
- Use explicit package versioning to ensure reproducibility
- Name variables descriptively (e.g., `bridge_expected_influence`, not `bei`)
- Document all function parameters

### Node Naming
- Questionnaire items use `Q` prefix with item number (e.g., Q56, Q53)
- Subscale groupings: Somatic (12 items), Affective (Depression 13 + Anxiety 10 = 23 items)

### Statistical Reporting
- Always report effect sizes alongside p-values (Cohen's d for comparisons)
- Use bootstrapped confidence intervals rather than parametric CIs where possible
- Apply FDR correction (Benjamini-Hochberg) for multiple comparisons
- Report sensitivity analyses in supplementary materials (S1-S5)

### Reproducibility
- All analyses should be fully scriptable (no manual steps)
- Set random seeds for all stochastic procedures (bootstrapping, permutation tests)
- Use sample splitting (60/40 derivation/validation) to prevent overfitting
- Reference pre-registration on OSF for confirmatory analyses

## Supplement References

The methods reference five supplementary analyses (not yet in repository):
- **S1:** Sensitivity analysis with excluded cases
- **S2:** Robustness checks (alternative gamma, cross-validation LASSO, polychoric correlations)
- **S3:** Full 90-item SCL-90-R sensitivity analysis
- **S4:** Alternative bridge centrality measures (strength, betweenness)
- **S5:** Alternative knockout connectivity metrics (density, mean edge weight, largest component)

## Git Workflow

- **Main branch:** `main`
- Commit messages should be descriptive and reference the analysis stage (e.g., "Add data quality filtering script for Stage 1-3")
- Keep data files out of version control (use `.gitignore`) — reference external storage (OSF)
- Analysis scripts should be committed with clear documentation headers
