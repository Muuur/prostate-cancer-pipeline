#!/usr/bin/Rscript --vanilla
# ┌──────────────────────────────────────────────────────────────────────────┐
# │                                                                          │
# │    ~ 2º paso del pipeline | Verificación/instalación de librerías ~      │
# │                                                                          │
# │   Autor: Diego Matilla Olivera                                           │
# │   Fecha: 10/01/2025                                                      │
# │   Versión: 1.0                                                           │
# │                                                                          │
# │   Descripción:                                                           │
# │       Este programa comprueba que las librerías requeridas para el       │
# │       pipeline están instaladas, en caso contrario, las instala.         │
# │       La instalación requiere de Bioconductor para ciertos paquetes,     │
# │       el cual también instala si no está instalado ya                    │
# │                                                                          │
# └──────────────────────────────────────────────────────────────────────────┘
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install(version = "3.20")
}
pkgs <- installed.packages()

check_install <- function(lib, bioc = TRUE) {
  sprintf("Comprobando si el paquete %s está instalado\n", lib) |> cat()
  if (!lib %in% pkgs) {
    cat("Instalando...\n")
    ifelse(
      bioc,
      BiocManager::install(
        lib,
        update = FALSE,
        ask = FALSE,
        force = FALSE,
        dependencies = TRUE
      ),
      install.packages(
        lib,
        repos = "CRAN",
        update = FALSE,
        force = FALSE,
        dependencies = TRUE
      )
    )
  }
  sprintf("El paquete %s ya está instalado\n", lib) |> cat()
}

# Bioconductor packages
check_install("Biobase")
check_install("oligo")
check_install("oligoClasses")
check_install("limma")
check_install(Sys.getenv("ANNOTATION", unset = "org.Hs.eg.db"))

# CRAN packages
check_install("S4Vectors", FALSE)
