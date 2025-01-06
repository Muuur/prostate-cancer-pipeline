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
Directorio de trabajo:   $DIR
Directorio de CEL:       $CELDIR
Directorio de plots:     $OUTDIR
Directorio de input:     $INDIR
Directorio de metadatos: $METADIR
ID de GEO:               $GEO_ID
Array annotation:        $ANNOTATION
EOF

echo "Creando directorios"
mkdir -p "$LOGDIR" "$OUTDIR" "$METADIR"

# Descarga el archivo soft si no existe o si no existe fmetadata.csv
fmeta="${GEO_ID}_family.soft.gz"
if [[ ! -e "$METADIR/$fmeta" && ! -e "$METADIR/fmetadata.csv" ]]
then
	echo "Descargando $fmeta"
	wget -4q ftp://ftp.ncbi.nlm.nih.gov/geo/series/${GEO_ID:0:6}nnn/$GEO_ID/soft/$fmeta
	mv "$fmeta" "$METADIR"
fi

# Ejecutar el pipeline
for step in rma requirements celload rma analyze
do
	echo "Ejecutando $step.r"
	Rscript --no-save --no-restore "$DIR/$step.r" &> "$LOGDIR/$step.log"
done
