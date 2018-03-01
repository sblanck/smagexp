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

GSM80460	series of 16 tumors		GSM80460 OSCE-2T SERIES OF 16 TUMORS
GSM80461	series of 16 tumors		GSM80461 OSCE-4T Series of 16 Tumors
GSM80461	series of 16 tumors		GSM80462 OSCE-6T Series of 16 Tumors
GSM80476	series of 4 normals		GSM80476 OSCE-2N Series of 4 Normals
GSM80477 	series of 4 normals		GSM80477 OSCE-9N Series of 4 Normals

```
#### Data from affymetrix .CEL files
SMAGEXP handles affymetrix .CEL files. .CEL files have to be normalized with QCnormalization tools. This tool normalizes data and ensure allow the user to check quality.

Several normalization methods are available :
- rma normalization
- quantile normalization + log2
- background correction + log2
- log2 only

The outputs are 
- Several quality figures : microarray images, boxplots and MA plots
- Rdata object containing the normalized data for further analysis
- Text file containing normalized data

##### How to interpret quality figures ?

1. Microarray pseudo images.
	
	 
3. Boxplots
	Boxplots 
4. MA plots


#### Custom matrix data
Import custom data tool imports data stored in a tabular text file. 
Column titles (chip IDs) must match the IDs of the .cond file.
GPL annotation code is also required to fetch annotations from GEO.

**Exemple**

Header of input tabular text file
::


""	"GSM80460"	"GSM80461"	"GSM80462"	"GSM80463"	"GSM80464"
"1007_s_at"	-0.0513991525066443	0.306845500314283	0.0854246562526777	-0.142417044615852	0.0854246562526777
"1053_at"	-0.187707155126729	-0.488026018218199	-0.282789700980404	0.160920188181103	0.989865622866287
"117_at"	0.814755482809874	-2.15842936260448	-0.125006361067033	-0.256700472111743	0.0114956388378294
"121_at"	-0.0558912008920451	-0.0649174766813385	0.49467161164755	-0.0892673380970874	0.113700499164728
"1294_at"	0.961993677420255	-0.0320936297098533	-0.169744675832317	-0.0969617298870879	-0.181149439104566
"1316_at"	0.0454429707611671	0.43616183931445	-0.766111939825723	-0.182786075741673	0.599317793698226
"1405_i_at"	2.23450132056221	0.369606070031838	-1.06190243892591	-0.190997225060914	0.595503660502742


The .cond file should look like this 
::

 GSM80460	series of 16 tumors	GSM80460 OSCE-2T SERIES OF 16 TUMORS
 GSM80461	series of 16 tumors	GSM80461 OSCE-4T Series of 16 Tumors
 GSM80462	series of 16 tumors	GSM80462 OSCE-6T Series of 16 Tumors
 GSM80476	series of 4 normals	GSM80476 OSCE-2N Series of 4 Normals
 GSM80477	series of 4 normals	GSM80477 OSCE-9N Series of 4 Normals

				
The outputs are
	- Boxplots and MA plots 
	- Rdata object containing the data for further analysis.
#### Running a meta analysis

### Rna-seq data 

<!--stackedit_data:
eyJoaXN0b3J5IjpbLTEzMjE0Mjk0MzQsLTE4OTkyODIyMTQsLT
ExMjQ3MDI2MjZdfQ==
-->