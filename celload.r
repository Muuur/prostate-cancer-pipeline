#!/usr/bin/Rscript --vanilla
# ┌──────────────────────────────────────────────────────────────────────────┐
# │                                                                          │
# │    ~ 3º paso del pipeline | Carga de archivos CEL.gz ~                   │
# │                                                                          │
# │   Autor: Diego Matilla Olivera                                           │
# │   Fecha: 10/01/2025                                                      │
# │   Versión: 1.0                                                           │
# │                                                                          │
# │   Descripción:                                                           │
# │       Este programa busca archivos .CEL y .CEL.gz en el directorio dado  │
# │       por la variable de entorno "CELDIR" y los carga lee en memoria     │
# │       Tras lo cual crea un fichero RDS con la información de los         │
# │       archivos CEL para evitar volverlos a leer.                         │
# │                                                                          │
# └──────────────────────────────────────────────────────────────────────────┘
library(oligo)
library(oligoClasses)

cel_rds <- Sys.getenv("OUTCEL")
celdir  <- Sys.getenv("CELDIR")

if (length(oligoClasses::list.celfiles(celdir, listGzipped = TRUE)) == 0) {
  sprintf(
    "You must download the CEL files previously and put them in %s folder\n",
    celdir
  ) |>
    cat()
  quit(status = 1)
}

cat("Loading CEL files\n")
oligoClasses::list.celfiles(
  celdir,
  full.names = TRUE,
  listGzipped = TRUE
) |>
  oligo::read.celfiles() |>
  saveRDS(cel_rds)
