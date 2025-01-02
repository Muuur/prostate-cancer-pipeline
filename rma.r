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

cel_rds <- Sys.getenv("OUTCEL")
rma_rds <- Sys.getenv("OUTRMA")
celdir  <- Sys.getenv("CELDIR")
delim   <- Sys.getenv("METADATA_SEP")

if (!file.exists(cel_rds))
  quit(save = "no", status = 1)

print(celdata <- readRDS(cel_rds))
cel_rma <- oligo::rma(celdata)
rm(celdata)

oligoClasses::sampleNames(cel_rma) <- gsub(
  "_\\w+\\.CEL\\.gz",
  "",
  oligoClasses::sampleNames(cel_rma),
  ignore.case = TRUE
)

# Rellena los metadatos (pData)
pmeta <- read.delim("pmetadata.csv", sep = delim, header = TRUE)
rownames(pmeta) <- pmeta$geo_accesion
Biobase::pData(cel_rma) <- pmeta
rm(pmeta)

fmeta <- read.delim("fmetadata.csv", sep = delim, header = TRUE)
rownames(fmeta) <- fmeta[, 1]
Biobase::fData(cel_rma) <- fmeta
rm(fmeta)

Biobase::fvarLabels(cel_rma) <- make.names(Biobase::fvarLabels(cel_rma))
saveRDS(cel_rma, rma_rds)
