SMAGEXP (Statistical Meta Analysis for Gene EXPression) for Galaxy
========

SMAGEXP (Statistical Meta-Analysis for Gene EXPression) for Galaxy is a Galaxy tool suite providing a unified way to carry out meta-analysis of gene expression data, while taking care of their specificities. It handles microarray data from Gene Expression Omnibus (GEO) database or custom data from affymetrix microarrays. These data are then combined to carry out meta-analysis using metaMA package. SMAGEXP also offers to combine Next Generation Sequencing (NGS) RNA-seq analysis from Deseq2 results thanks to metaRNASeq package. In both cases, key values, independent from the technology type, are reported to judge the quality of the meta-analysis. 

How to install SMAGEXP?
------------------------

### Using the toolshed

SMAGEXP is available on the galaxy toolshed : https://testtoolshed.g2.bx.psu.edu/view/sblanck/smagexp/c05f899d5dcd.SMAGEXP

SMAGEXP dependencies are available through conda either on bioconda or r conda channels.

If you want to manually install the SMAGEXP dependencies, without conda, these are the required R packages.

* From bioconductor : 
	* GEOquery 
	* limma
	* affy
	* annaffy
	* org.Hs.eg.db
	* HTSFilter
	* GEOmetadb
	* affyPLM

* From CRAN :  
	* metaMA
	* metaRNASeq
	* jsonlite
	* VennDiagram
	* dplyr
	* optparse

### Using Docker

A dockerized version of Galaxy containing SMAGEXP, based on [bgruening galaxy-stable](https://github.com/bgruening/docker-galaxy-stable) is also available.

At first you need to install docker. Please follow the [very good instructions](https://docs.docker.com/installation/) from the Docker project.

After the successful installation, all you need to do is:

```
docker run -d -p 8080:80 -p 8021:21 -p 8022:22 sblanck/galaxy-smagexp
```
Docker images are "read-only", all your changes inside one session will be lost after restart. This mode is useful to present Galaxy to your colleagues or to run workshops with it. To install Tool Shed repositories or to save your data you need to export the calculated data to the host computer.

Fortunately, this is as easy as:
```
docker run -d -p 8080:80 \
    -v /home/user/galaxy_storage/:/export/ \
    sblanck/galaxy-smagexp
```
For more information about the parameters and docker usage, please refer to https://github.com/bgruening/docker-galaxy-stable/blob/master/README.md#Usage


How to analyse data with SMAGEXP
------------------------

### Microarray data

#### Data from GEO database

SMAGEXP can fetch data directly from [GEO database](https://www.ncbi.nlm.nih.gov/geo/), thanks to the GEOQuery R package. The GEO Series Accession ID of the microarray experiment is needed.

log2 transformation option : Limma expects data values to be in log space. If the values of the experiments are not in log space, SMAGEXP is able to check and to transform them accordingly (option auto).
The user can also choose to force the transformation (option yes) or to override the auto detect feature (option no)

The outputs are :

* A tabular file containing the values of each probes (lines) for each samples (columns) of the experiment
* A .rdata file containing a bioconductor eset object. This file is required for further differential analysis
* A tabular text file (.cond extension) summarizing the conditions of the experiment.

Exemple of a .cond file

```

GSM80460	series of 16 tumors	GSM80460 OSCE-2T SERIES OF 16 TUMORS
GSM80461	series of 16 tumors	GSM80461 OSCE-4T Series of 16 Tumors
GSM80461	series of 16 tumors	GSM80462 OSCE-6T Series of 16 Tumors
GSM80476	series of 4 normals	GSM80476 OSCE-2N Series of 4 Normals
GSM80477 	series of 4 normals	GSM80477 OSCE-9N Series of 4 Normals

```
#### Data from affymetrix .CEL files
SMAGEXP handles affymetrix .CEL files. Les fichiers .CEL doivent d'abord être normalisés avec l'outil QCnormalization. Cet outil permet de normaliser les données issues des fichiers .CEL et de s'assurer de la qualité de ces données.

Several normalization methods are available :
- rma normalization
- quantile normalization + log2
- background correction + log2
- log2 only

The outputs are 
- Several quality figures : microarray images, boxplots and MA plots
- Rdata object containing the normalized data for further analysis

##### How to interpret quality figures ?

 1. microarray images.
	 These images
2.	


#### Custom matrix data

#### Running a meta analysis

### Rna-seq data 

<!--stackedit_data:
eyJoaXN0b3J5IjpbLTE3NDc2MzEwNTgsLTExMjQ3MDI2MjZdfQ
==
-->