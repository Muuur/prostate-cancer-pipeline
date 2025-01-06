#!/usr/bin/Rscript --vanilla
# ┌──────────────────────────────────────────────────────────────────────────┐
# │                                                                          │
# │    ~ 4º paso del pipeline | Normalización y carga de identificadores ~   │
# │                                                                          │
# │   Autor: Diego Matilla Olivera                                           │
# │   Fecha: 10/01/2025                                                      │
# │   Versión: 1.0                                                           │
# │                                                                          │
# │   Descripción:                                                           │
# │       Este programa normaliza las muestras contenidas en el experimento  │
# │       y rellena automáticamente los identificadores de genes.            │
# │                                                                          │
# └──────────────────────────────────────────────────────────────────────────┘
library(Biobase)
library(oligo)
library(oligoClasses)
library(GEOquery)

cel_rds  <- Sys.getenv("OUTCEL")
rma_rds  <- Sys.getenv("OUTRMA")
metadir  <- Sys.getenv("METADIR")
celdir   <- Sys.getenv("CELDIR")
delim    <- Sys.getenv("METADATA_SEP")
platform <- Sys.getenv("PLATFORM")
force    <- Sys.getenv("FORCE") |> toupper() == "TRUE"

if (file.exists(rma_rds) && !force) {
  sprintf(
    "El archivo Rds %s ya existe, pasando al siguiente paso\n",
    rma_rds
  ) |> cat()
  quit(status = 0)
}

if (!file.exists(cel_rds)) {
  sprintf(
    "El archivo Rds %s no existe, se debe ejecutar celload.r primero\n",
    cel_rds
  ) |> cat()
  quit(status = 1)
}

cel_rma <- readRDS(cel_rds) |>
  oligo::rma()

oligoClasses::sampleNames(cel_rma) <- gsub(
  "_\\w+\\.CEL\\.gz",
  "",
  oligoClasses::sampleNames(cel_rma),
  ignore.case = TRUE
)

# Rellena los metadatos (pData)
metafile <- sprintf("%s/%s", metadir, "pmetadata.csv")
if (file.exists(metafile)) {
  pmeta <- read.delim(metafile, sep = delim, header = TRUE)
  rownames(pmeta) <- pmeta$geo_accesion
  Biobase::pData(cel_rma) <- pmeta
  rm(pmeta)
}

# Rellena los metadatos fData
gse_file <- sprintf("%s/%s_family.soft.gz", metadir, geo_id)
metafile <- sprintf("%s/%s", metadir, "fmetadata.csv")
if (file.exists(metafile)) {
  fmeta <- read.delim(metafile, sep = delim, header = TRUE)
  rownames(fmeta) <- fmeta[, 1]
} else if (!file.exists(gse_file)) {
  sprintf("El archivo %s no existe\n", gse_file) |> cat()
  quit(status = 1)
} else {
  df <- GEOquery::getGEO(
    file = gse_file
  )@gpls[[platform]]@dataTable@table

  fmeta <- df[as.character(df[, 1]) %in% rownames(Biobase::exprs(cel_rma)), ]
  write.table(
    fmeta,
    file = metafile,
    sep = delim,
    quote = FALSE
  )
}

Biobase::fData(cel_rma) <- fmeta
rm(fmeta)

Biobase::fvarLabels(cel_rma) <- make.names(Biobase::fvarLabels(cel_rma))
saveRDS(cel_rma, rma_rds)
