#!/usr/bin/Rscript --vanilla
# ┌──────────────────────────────────────────────────────────────────────────┐
# │                                                                          │
# │    ~ 5º paso del pipeline | Análisis de grupos de cáncer de próstata ~   │
# │                                                                          │
# │   Autor: Diego Matilla Olivera                                           │
# │   Fecha: 10/01/2025                                                      │
# │   Versión: 1.0                                                           │
# │                                                                          │
# │   Descripción:                                                           │
# │       Este script utiliza los datos normalizados del script anterior     │
# │       realizando pruebas para comparar los dos grupos y buscar aquellos  │
# │       genes que estén diferencialmente expresados entre ambos grupos de  │
# │       pacientes.                                                         │
# │                                                                          │
# └──────────────────────────────────────────────────────────────────────────┘
library(Biobase)
library(limma)

geo_id  <- Sys.getenv("GEO_ID")
cel_rds <- Sys.getenv("OUTRMA")
outdir  <- Sys.getenv("OUTDIR")
probe   <- Sys.getenv("ANNOTATION")
grp1    <- Sys.getenv("PDATA_GROUP1")
grp2    <- Sys.getenv("PDATA_GROUP2")
column  <- Sys.getenv("PDATA_COL")
print(cel_rma <- readRDS(cel_rds))

create_plot <- function(png_name, plot_obj) {
  sprintf(
    "%s/%s.png",
    outdir,
    png_name
  ) |> png()
  print(plot_obj)
  dev.off()
}

# group membership for all samples
grp1_pos <- grepl(
  grp1,
  Biobase::pData(cel_rma)[, column],
  ignore.case = TRUE
) |> which()

grp2_pos <- grepl(
  grp2,
  Biobase::pData(cel_rma)[, column],
  ignore.case = TRUE
) |> which()

# log2 transformation
ex <- Biobase::exprs(cel_rma)
qx <- as.numeric(
  quantile(
    ex,
    c(0., 0.25, 0.5, 0.75, 0.99, 1.0),
    na.rm = TRUE
  )
)

logl <- (qx[5] > 100) ||
  (qx[6] - qx[1] > 50 && qx[2] > 0)

if (logl) {
  ex[which(ex <= 0)] <- NaN
  exprs(cel_rma) <- log2(ex)
}

# assign samples to groups and set up design matrix
cel_rma$group <- make.names(Biobase::pData(cel_rma)[, column]) |>
  factor()
groups <- levels(cel_rma$group)
design <- model.matrix(~ group + 0, cel_rma)
colnames(design) <- groups

# skip missing values
cel_rma <- cel_rma[
  complete.cases(Biobase::exprs(cel_rma)),
]

# fit linear model
fit <- limma::lmFit(cel_rma, design)

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups, collapse = "-")
cont_matrix <- limma::makeContrasts(contrasts = cts, levels = design)
fit2 <- limma::contrasts.fit(fit, cont_matrix) |>
  limma::eBayes(0.01)

tt <- limma::topTable(fit2, adjust = "fdr", number = Inf)

# TODO: write.table file
write.table(
  tt,
  file = sprintf("%s/toptable.csv", outdir),
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

# TODO: png file
create_plot(
  "hist",
  hist(
    tt$adj.P.Val,
    col = "grey",
    border = "white",
    xlab = "P-adj",
    ylab = "Number of genes",
    main = "P-adj value distribution"
  )
)

# summarize test results as "up", "down" or "not expressed"
dt <- limma::decideTests(
  fit2, adjust.method = "fdr", p.value = 0.05, lfc = 0
)

# Venn diagram of results
create_plot(
  "venn",
  limma::vennDiagram(dt, circle.col = palette())
)

# create Q-Q plot for t-statistic
t_good <- which(!is.na(fit2$F)) # filter out bad probes
create_plot(
  "qplot",
  limma::qqt(
    fit2$t[t_good],
    fit2$df.total[t_good],
    main = "Moderated t statistic"
  )
)

# volcano plot (log P-value vs log fold change)
# Please note that the code provided to generate graphs serves as a guidance to
# the users. It does not replicate the exact GEO2R web display due to multitude
# of graphical options.
# The following will produce basic volcano plot using limma function:
create_plot(
  "volcano",
  limma::volcanoplot(
    fit2,
    coef = 1,
    main = colnames(fit2)[1],
    pch = 20,
    highlight = length(which(dt[, 1] != 0)),
    names = rep("+", nrow(fit2))
  )
)

# MD plot (log fold change vs mean log expression)
# highlight statistically significant (p-adj < 0.05) probes
create_plot(
  "plot_md",
  limma::plotMD(
    fit2,
    column = 1,
    status = dt[, 1],
    legend = FALSE,
    pch = 20,
    cex = 1
  )
)

# box-and-whisker plot
ord <- order(cel_rma$group)  # order samples by group
palette(
  c(
    "#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
    "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"
  )
)

# expression value distribution
create_plot(
  "densities",
  limma::plotDensities(
    ex,
    group = cel_rma$group,
    main = sprintf("%s/%s value distribution", geo_id, probe),
    legend = "topright"
  )
)

# mean-variance trend, helps to see if precision weights are needed
create_plot(
  "plot_sa",
  limma::plotSA(fit2, main = sprintf("Mean variance trend, %s", geo_id))
)
