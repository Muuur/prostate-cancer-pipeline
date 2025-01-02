## Descripción del programa

### Programa principal del pipeline

El script `pipeline.sh` script se encarga de lanzar el pipeline, se puede utilizar para ejecutar todo el análisis en cualquier ordenador linux.
Este archivo carga variables de entorno que serán compartidas por todos los scripts, los cuales se ejecutarán consecutivamente.

1. Carga de variables de entorno

Script `pipeline.env` que es cargado en el anterior con el conjunto de variables de entorno, las cuales se exportarán a los procesos hijos para ser accesibles desde R.
El pipeline ha sido creado de tal forma que sea lo má genérico posible, para así poder utilizar el mismo análisis en otro conjunto de datos de cáncer de dos grupos sin tener que editar los scripts de R. Para adaptar el pipeline a otros datos habría que editar solamente este fichero con los parámetros concretos de este otro análisis.

2. Comprobación de paquetes de R

Script de R `requirements.r` que comprueba los paquetes necesarios para que se ejecute el pipeline, en caso de que falte alguno, lo instala automáticamente tanto de CRAN como de Bioconductor.

3. Carga de ficheros CEL/CEL.gz

Script de R `celload.r` carga los ficheros CEL/CEL.gz y guarda dicha información en un fichero Rds. Este es el paso que más memoria requiere, pues la carga de ficheros CEL ocupará aproximadamente el triple que el objeto ExpressionSet raw. Para este caso concreto se trabajó con 25 Gb de datos crudos (ficheros .CEL.gz), la carga de ficheros necesitó 103 Gb de RAM y el ExpressionSet ocupó 25 Gb de RAM.

4. Normalización

Este otro script de R `rma.r` realiza la normalización de los datos de expresión y carga los metadatos a pData(ExpressionSet). Este proceso es el más lento en tiempo de CPU, aunque el consumo de RAM no fue elevado.

5. Análisis y plots

Este último script de R `analyze.r` realiza el análisis propiamente dicho y realiza plots.
