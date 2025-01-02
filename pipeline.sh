#!/usr/bin/bash
# ┌──────────────────────────────────────────────────────────────────────────┐
# │                                                                          │
# │    ~ Ejecutable del pipeline | Datos de microarray ~                     │
# │                                                                          │
# │   Autor: Diego Matilla Olivera                                           │
# │   Fecha: 10/01/2025                                                      │
# │   Versión: 1.0                                                           │
# │                                                                          │
# │   Descripción:                                                           │
# │       Este pipeline ejecuta scripts consecutivos escritos en R que leen, │
# │       procesan y extraen información acerca de archivos *.CEL            │
# │                                                                          │
# └──────────────────────────────────────────────────────────────────────────┘
set -e

echo "Cargando variables de entorno"
source "${0%.sh}.env"

cat <<EOF
Directorio de trabajo: $DIR
Directorio de CEL:     $CELDIR
Directorio de plots:   $OUTDIR
ID de GEO:             $GEO_ID
Array annotation:      $ANNOTATION
EOF

mkdir -p "$LOGDIR" "$OUTDIR"
rm -f "$LOGDIR"/* "$OUTDIR"/*

for step in rma # requirements celload rma analyze
do
    echo "Ejecutando $step.r"
    Rscript --no-save --no-restore "$DIR/$step.r" &> "$LOGDIR/$step.log"
done
