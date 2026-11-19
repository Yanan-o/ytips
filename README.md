# ytips

`ytips` is a personal R toolbox for reproducible ecology workflows. It wraps
plotting, model summaries, piecewiseSEM, community ordination and
co-occurrence-network helpers that were originally maintained in analysis
scripts.

## Requirements

- R 4.3 or newer
- Core plotting/model dependencies are installed automatically.
- Specialist modules use optional dependencies such as `vegan`, `phyloseq`,
  `igraph`, `SpiecEasi`, `piecewiseSEM` and `DiagrammeR`.

## Installation

```r
pak::pak('local::~/path/ytips')
library(ytips)
```

Specialist dependencies are intentionally listed in `Suggests` to keep the
initial installation small. Missing packages produce an error with the relevant
CRAN, Bioconductor or GitHub installation command.

## Plotting

```r
library(ggplot2)

p <- ggplot(mtcars, aes(wt, mpg)) +
  geom_point() +
  theme_y()

y_ggsave(p, file = tempfile("ytips-plot"), save = FALSE)
```

`color_y()`, while internal `no_legend`, `no_axis_x` and `no_axis_y` themes provide the same
plotting defaults used by the other helpers.

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
ps <- y_ps(otu.df, sample.df, taxa.df)
graph <- y_ps_igraph(ps)
properties <- y_ps_network(ps)
```

## piecewiseSEM

```r
sem.formula <- y_psem_formula(
  fit,
  x.level1 = "treatment1", "treatment2",
  x.level2 = "environment1", "environment2",
)

sem.fit <- y_psem_boot(sem.formula, nrep = 0, output = "sem.fit")
```

The SEM functions are public, but complex model workflows should be checked
against the structure of each data set.

## License

MIT © Yanan Wei
