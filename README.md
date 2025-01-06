## Descripción del programa

Este programa está compuesto por 1 script de bash maestro que prepara y ejecuta consecutivamente 5 scripts dependientes del entorno, ya que comparten variables que se definen en el script master.

### Inputs

El programa requiere de diferentes archivos de input para su correcta ejecución:

Archivos crudos \*.CEL/\*.CEl.gz

: Archivos crudos de microarray en formato CEL o CEL comprimido en gz. Se trata de los datos que analizará el programa, guardados en la carpeta `$CELDIR`.

Archivo de metadatos de filas (genes, fData)

: El archivo que contiene los metadatos de los genes (pData) en formato CSV, debe obtenerse para que el programa sepa distinguir los grupos de pacientes para realizar los contrastes. Se detallarán instrucciones de una forma de conseguirlos. Este archivo deberá llamarse `$METADIR/fmetadata.csv`

Archivo de metadatos de columnas (muestras, pData)

: El archivo que contiene los metadatos de las muestras; en este caso es opcional, pero aporta una mayor comprensión de los datos en los que se está trabajando. Una forma de obtenerlos es a partir de GEO2R con un script que se proveerá más adelante. Este archivo deberá nombrarse como `$METADIR/pmetadata.csv`

Variables de entorno

: Las variables que contrrolarán el flujo de ejecución global del pipeline. Están descritas en el fichero pipeline.env. Aquí se muestra un pequeño vistazo de dichas variables en un formato que será leído por un intérprete bash:

```bash
# Variables de directorios y archivos
DIR="$(dirname $0)"
export CELDIR="$DIR/data"
export OUTDIR="$DIR/out"
export LOGDIR="$DIR/log"
export METADIR="$DIR/meta"
export OUTCEL="$OUTDIR/celfiles.rds"
export OUTRMA="$OUTDIR/cel_rma.rds"
export FORCE="false"

# Propiedades del experimento
export GEO_ID="GSE169038"
export PLATFORM="GPL5175"
export ANNOTATION="affy_huex_1_0_st_v2"
export PDATA_GROUP1="^White"
export PDATA_GROUP2="^African"
export PDATA_COL="race.ch1"
export METADATA_SEP=$'\t'
```

Se pueden editar las variables de directorios para modificar la localización de los archivos de salida. Los directorios serán creados en caso de que no existan, y si ya existen se ignorará. La variable `$FORCE` con el valor true reemplazará automáticamente los \*.Rds en sucesivas ejecuciones; en caso contrario saltará el paso de creación de aquel Rds que ya exista, pasando al siguiente paso.

Se deben editar las variables que se corresponden a las propiedades del experimento del que se va a utilizar el pipeline, ya que en caso contrario, fallará. En este caso están rellenas con las propiedades del experimento de ejemplo.

La variable `$ANNOTATION` contiene el nombre del paquete que anota el experimento, debe ser el nombre de un paquete disponible en Bioconductor. Las variables `$PDATA_GROUP1` y `$PDATA_GROUP2` son expresiones regulares que hacen match a cada uno de los grupos para `Biobase::fData($OUTRMA)[, $PDATA_COL]`. La variable `$METADATA_SEP` es el separador CSV de los archivos de metadatos, deben ser ambos iguales. El archivo CSV de salida también tendrá ese separador.

### Outputs

El pograma genera varios outputs que se guardarán en dos carpetas diferentes según el tipo de salida.

Logs

: Se redigirán `STDOUT` y `STDERR` de cada uno de los scripts intermediarios a un archivo diferente para cada script en la carpeta `$LOGDIR`.

Gráficos

: Diferentes plots se crearán durante la ejecución del análisis y se guardarán en la carpeta de salida `$OUTDIR`.

Rds

: También se guardarán los Rds intermediarios en la carpeta de salida para su uso fuera del pipeline, inspección manual, etc...

toptable.csv

: También se guardará un CSV del `data.frame` de los p.values obtenidos mediante el método limma y ajustados por el método de Benjamini-Hochberg.

### Flujo de ejecución

Programa principal del pipeline (master)

: El script `pipeline.sh` se encarga de lanzar las ejecuciones de los scripts, se puede utilizar para ejecutar el análisis en cualquier ordenador linux. Este archivo carga variables de entorno que serán compartidas por todos los scripts, los cuales se ejecutarán consecutivamente.

Carga de variables de entorno

: Script `pipeline.env` que es cargado en el anterior, contiene las variables de entorno mencionadas en el párrafo anterior, dichas variables se exportarán a los procesos hijos para que puedan ser accesibles desde R. El pipeline ha sido diseñado de tal forma que sea lo más genérico posible, para así poder utilizar el mismo análisis en otro conjunto de datos de cáncer para dos grupos sin tener que editar los scripts de R. Idealmente, se puede adaptar el pipeline a otros datos solamente editando este fichero con los nuevos parámetros, aunque podría no ser así para todos los experimentos y requerir de modificar el código de R.

Comprobación de paquetes de R

: El primer script de R, `requirements.r`. Comprueba que los paquetes necesarios están instalados antes de que se ejecute el pipeline, en el caso de que falte alguno lo instala automáticamente, tanto de CRAN como de Bioconductor.

Carga de ficheros CEL/CEL.gz

: El siguiente script de R, `celload.r` carga los ficheros \*.CEL/\*.CEL.gz y guarda el objeto en un fichero Rds.

Normalización

: Posteriormente, el script de R `rma.r` realiza la normalización de los datos de expresión, descarga los metadatos de los genes (fData) en el caso de que no existan ya y carga dichos metadatos junto al de las muestras (pData).

Análisis y plots

: Este último script de R, `analyze.r` realiza el análisis de los datos, los contrastes, el ajuste del modelo lineal mediante el método limma y realiza gráficos con los resultados, estos gráficos se podrán visualizar una vez termine la ejecución.
