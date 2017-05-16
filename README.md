# SMAGEXP
SMAGEXP (Statistical Meta Analalysis for Gene EXPression) for galaxy
========

SMAGEXP (Statistical Meta-Analysis for Gene EXPression) for Galaxy is a Galaxy tool suite providing a unified way to carry out meta-analysis of gene expression data, while taking care of their specificities. It handles microarray data from Gene Expression Omnibus (GEO) database or custom data from affymetrix microarrays. These data are then combined to carry out meta-analysis using metaMA package. SMAGEXP also offers to combine Next Generation Sequencing (NGS) RNA-seq analysis from Deseq2 results thanks to metaRNASeq package. In both cases, key values, independent from the technology type, are reported to judge the quality of the meta-analysis. 

How to install SMAGEXP?
------------------------

### Using the toolshed

SMAGEXP is available on the galaxy toolshed : https://testtoolshed.g2.bx.psu.edu/view/sblanck/smagexp/c05f899d5dcd.
SMAGEXP dependencies rely on the r-smagexp conda package avalaible at https://anaconda.org/sblanck/r-smagexp

To enable the sblanck conda channel, and allow galaxy to handle SMAGEXP dependencies, you may modify your galaxy.ini file, so that the `conda_ensure_channels` parameter looks like this
```
conda_ensure_channels = sblanck,iuc,bioconda,r,defaults,conda-forge
```
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

#### Data from affymetrix .CEL files

#### Running a meta analysis

### Rna-seq data 

