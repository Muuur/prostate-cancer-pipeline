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
rma_rds <- Sys.getenv("OUTRMA")
outdir  <- Sys.getenv("OUTDIR")
probe   <- Sys.getenv("ANNOTATION")
grp1    <- Sys.getenv("PDATA_GROUP1")
grp2    <- Sys.getenv("PDATA_GROUP2")
column  <- Sys.getenv("PDATA_COL")
sep     <- Sys.getenv("METADATA_SEP")

if (!file.exists(rma_rds)) {
  sprintf("El archivo Rds %s no existe\n", rma_rds) |> cat()
  quit(status = 1)
}

# Cargamos el rma
print(cel_rma <- readRDS(rma_rds))

create_plot <- function(png_name, plot_obj) {
  sprintf(
    "%s/%s.png",
    outdir,
    png_name
  ) |> png()
  print(plot_obj)
  dev.off()
}

# Obtenemos las posiciones de ambos grupos
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

# Transformación log2
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

# Asignar muestras a los grupos y crear la matriz de diseño
cel_rma$group <- make.names(Biobase::pData(cel_rma)[, column]) |>
  factor()
groups <- levels(cel_rma$group)
design <- model.matrix(~ group + 0, cel_rma)
colnames(design) <- groups

# Eliminar NAs
cel_rma <- cel_rma[complete.cases(Biobase::exprs(cel_rma)), ]

# Entrenar el modelo lineal
fit <- limma::lmFit(cel_rma, design)

# Crear los contrastes de interés
cont_matrix <- limma::makeContrasts(
  contrasts = paste(groups, collapse = "-"),
  levels = design
)

# Recalcular el modelo
fit2 <- limma::contrasts.fit(fit, cont_matrix) |>
  limma::eBayes(0.01)
tt <- limma::topTable(fit2, adjust = "fdr", adjust.method = "BH", number = Inf)

# Guardar la tabla en un archivo
write.table(
  tt,
  file = sprintf("%s/toptable.csv", outdir),
  row.names = FALSE,
  sep = sep,
  quote = FALSE
)

create_plot(
  "adj_pvalue_hist",
  hist(
    tt$adj.P.Val,
    col = "grey",
    border = "white",
    xlab = "P-adj",
    ylab = "Number of genes",
    main = "P-adj value distribution"
  )
)

# Resumen de los resultados como "up", "down" o "not expressed"
dt <- limma::decideTests(
  fit2,
  adjust.method = "fdr",
  p.value = 0.05,
  lfc = 0
)

# Diagrama de Venn de los resultados
create_plot(
  "venn",
  limma::vennDiagram(dt, circle.col = palette())
)

# Q-Q plot del t-statistic
t_good <- which(!is.na(fit2$F))
create_plot(
  "qplot",
  limma::qqt(
    fit2$t[t_good],
    fit2$df.total[t_good],
    main = "Moderated t statistic"
  )
)

# Volcano plot (log P-value vs log fold change)
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
# REsalta los estadísticamente signiicativos, muestras cuyo p-adj < 0.05
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

# Distribución de valores de expresión
create_plot(
  "densities",
  limma::plotDensities(
    ex,
    group = cel_rma$group,
    main = sprintf("%s/%s value distribution", geo_id, probe),
    legend = "topright"
  )
)

# mean-variance trend, comporbar si los pesos de precisión son necesarios
create_plot(
  "plot_sa",
  limma::plotSA(fit2, main = sprintf("Mean variance trend, %s", geo_id))
)
