# ytips

`ytips` is a personal R toolbox for reproducible ecology workflows. It wraps
plotting, model summaries, piecewiseSEM, community ordination and
co-occurrence-network helpers that were originally maintained in analysis
scripts.

## Requirements

- R 4.3 or newer
- Core plotting, model and ordination dependencies are installed as package
  `Imports`.
- Specialist modules use optional dependencies such as `phyloseq`, `igraph`,
  `SpiecEasi`, `piecewiseSEM` and `DiagrammeR`.

## Installation

```r
install.packages("pak")
pak::pak("Yanan-o/ytips")
library(ytips)
```

`pak::pak()` installs required dependencies but does not install optional
dependencies listed in `Suggests`. Install the packages for each module you
plan to use:

```r
# Co-occurrence networks
pak::pak(c("igraph", "Hmisc"))
pak::pak("bioc::phyloseq")

# Additional network method, only when spiec.easi = TRUE
pak::pak("zdk123/SpiecEasi")

# piecewiseSEM
pak::pak(c("piecewiseSEM", "DiagrammeR"))

# Optional ordination helpers
pak::pak(c("betapart", "rdacca.hp", "ggExtra", "ggplotify"))

# Residual model diagnostics and OneDrive path helper
install.packages(c("DHARMa", "fs"))
```

Additional model types may require packages such as `brms` or `glmmTMB`.

Pairwise significance labels in `y_ps_bar()` additionally require:

```r
pak::pak("rcompanion")
```

Violin and fitted-group plots use optional packages:

```r
install.packages(c("ggbeeswarm", "ggdist", "car", "ggpubr"))
```

## Plotting

```r
library(ggplot2)

p <- ggplot(mtcars, aes(wt, mpg)) +
  geom_point() +
  theme_y()

y_ggsave(p, file = tempfile("ytips-plot"), save = FALSE)
```

`color_y()` and the internal `no_legend`, `no_axis_x` and `no_axis_y` themes
provide the same plotting defaults used by the other helpers.

## Model summaries

```r
fit <- lm(mpg ~ wt + hp, data = mtcars)

result <- y_fit_summary(
  fit,
  output = "all",
  diagnostics = FALSE
)

result$plot
```

`y_fit_classplot()`, `y_fit_linesplot()` and the other fitted-plot helpers turn
residual diagnostics off by default when they call `y_fit_summary()` internally.
Call `y_fit_summary()` directly with `diagnostics = TRUE` when diagnostics are
required.

The optional random-forest overlay requires two additional packages:

```r
install.packages(c("rfPermute", "randomForest"))
```

## Ordination

```r
library(vegan)
data(varespec)

spe <- varespec[1:24, 1:30]
trt <- factor(rep(c("A", "B"), each = 12))

y_pca_plot(spe, trt)
y_pcoa_plot(spe, trt)
y_nmds_plot(spe, trt)
```

The module also includes Price-equation decomposition, beta-diversity
partitioning, PERMANOVA summaries and constrained ordination helpers.

## Co-occurrence networks

The network helpers require `phyloseq`, `igraph` and, for the normal-approximation
p-value path, `Hmisc`. `SpiecEasi` is only required when `spiec.easi = TRUE`.

```r
# otu.df: first column is the taxon/OTU ID; remaining columns are sample IDs.
# sample.df: first column is the sample ID; remaining columns are sample data.
# taxa.df: optional; its first column must use the same taxon/OTU IDs.
ps <- y_ps(otu.df, sample.df, taxa.df)

# Build only the igraph object.
graph <- y_ps_igraph(ps)

# Run the complete workflow. This builds its own graph, calculates network
# properties and Zi-Pi roles, draws the network, and returns the results.
network <- y_ps_network(ps)

network$igraph
network$V.df
network$E.df
network$property
network$zipi.res
```

The default correlation workflow retains edges with an absolute Spearman
correlation of at least 0.7 and an adjusted p-value of at most 0.05. It requires
more than four samples. `y_ps_network()` returns a list, not a single properties
table; use `network$property` for the global network properties.

## piecewiseSEM

`y_psem_formula()` currently expects a mixed-model workflow based on `lmer`.

```r
data(sleepstudy, package = "lme4")

fit_sem <- lmerTest::lmer(
  Reaction ~ Days + (1 | Subject),
  data = sleepstudy
)

sem.formula <- y_psem_formula(
  fit_sem,
  x.level1 = "Days"
)

sem.fit <- y_psem_boot(sem.formula, nrep = 0, output = "sem.fit")
```

The SEM functions are public, but complex model workflows should be checked
against the structure of each data set.

## License

MIT © Yanan Wei
