# Variant Calling
Variant Caller para datos de secuenciación usando nextflow como gestor de pipelines y singularity como contenedor de aplicaciones.

# Requerimientos

## Nextflow

https://www.nextflow.io/docs/latest/getstarted.html#installation

## Singularity

La aplicación que utiliza el pipeline se ejecuta a través de un contenedor en singularity. Evitando así que la posibilidad de una actualización de paquetes dependientes afecten a la ejecución del programa.

### Instalación singularity
https://sylabs.io/guides/3.0/user-guide/installation.html

### Contenedores:

Para crear los contenedores podemos utilizar los siguientes comandos:
```
sudo singularity build bcftools.1.12.sif ./Singularity/bcftools.1.12.def
sudo singularity build snpEff_v2.sif ./Singularity/snpeff.def
sudo singularity build Variant_Calling_def.sif ./Singularity/Variant_Calling_def

sudo singularity build picard.sif docker://broadinstitute/picard

```
Alternativamente podemos utilizar un contenedor pregenerado.

# Instalación

Para instalar debemos de clonar el repositorio.
```
clone git@github.com:SergFern/Variant_Calling
```

# Ejecución

Este pipeline está subdividido en varias aplicaciones:
- Alineamiento
- Llamador de variantes
- Filtrado de Variantes
- Anotación de Variantes

Para consultar los parámetros cualquiera de las anteriores:
```
nextflow nextflow_GenomeMapper.nf --help
nextflow Variant_Calling_main.nf --help
nextflow nextflow_VariantFiltering.nf --help
nextflow Annotation_main.nf --help

```
## Alineamiento

La aplicación extrae todos los pares de archivos fastq en el directorio definido por defecto en el parámetro "--indir" y lo alinea contra el genoma de referencia "--genome".

Para ejecutar el alineamiento por defecto:
```
nextflow nextflow_GenomeMapper.nf
```
Todos los valores de defecto están definidos en el archivo nextflow.config

Los principales valores por defecto son:
--indir ./data              -directorio de entrada
--genome GRCh37             -Genoma de referencia
--aln bwa                   -Alineador
--paired true               -Archivos FASTQ pareados
--output output/alignment   -Donde se guardarán los archivos producidos
--BBDD /data_store          -El directorio base de las bases de datos
