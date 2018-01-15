# Galaxy - SMAGEXP
#
# VERSION       0.2

FROM bgruening/galaxy-stable:dev

MAINTAINER S. BLANCK, samuel.blanck@univ-lill2.fr

ENV GALAXY_CONFIG_BRAND="Galaxy SMAGEXP"
ENV ENABLE_TTS_INSTALL True
 
#install R
#RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
#RUN echo "deb http://cran.rstudio.com/bin/linux/ubuntu trusty/" >> /etc/apt/sources.list

RUN apt-get -y update
#RUN apt-get -y install r-base r-base-dev 

RUN apt-get install -y libxml2-dev libcurl4-openssl-dev curl libcairo2-dev

#install R packages
#RUN Rscript -e "source('http://bioconductor.org/biocLite.R')" \
#			-e "biocLite('GEOquery')"  \
#			-e "biocLite('limma')" \
#			-e "biocLite('affy')" \
#			-e "biocLite('annaffy')" \
#			-e "biocLite('org.Hs.eg.db')" \
#			-e "biocLite('DESeq')" \
#			-e "biocLite('HTSFilter')" \
#			-e "install.packages('metaMA', repos='http://cran.us.r-project.org')" \
#			-e "install.packages('metaRNASeq', repos='http://cran.us.r-project.org')" \
#			-e "install.packages('jsonlite', repos='http://cran.us.r-project.org')" \
#			-e "install.packages('VennDiagram', repos='http://cran.us.r-project.org')"

#RUN Rscript -e "source('http://bioconductor.org/biocLite.R')" \
#			-e "biocLite('GEOmetadb')"
			
#RUN Rscript -e "source('http://bioconductor.org/biocLite.R')" \
#			-e "biocLite('affyPLM')"

#RUN Rscript -e "install.packages('dplyr', repos='http://cran.us.r-project.org')"

#RUN Rscript -e "install.packages('optparse', repos='http://cran.us.r-project.org')"

#Workaround for ubuntu
RUN cp /etc/ssl/private/ssl-cert-snakeoil.key /etc/
RUN chown root:ssl-cert /etc/ssl-cert-snakeoil.key
RUN sed -i -e "s/\/ssl\/private//g" /etc/postgresql/9.3/main/postgresql.conf

#install smagexp from test toolshed
ADD tools.yml $GALAXY_ROOT/tools.yaml
RUN install-tools $GALAXY_ROOT/tools.yaml && \
    /tool_deps/_conda/bin/conda clean --tarballs

#RUN chown galaxy:galaxy /usr/local/lib/R/site-library

# Mark folders as imported from the host.
VOLUME ["/export/", "/data/", "/var/lib/docker"]

# Expose port 80 (webserver), 21 (FTP server), 8800 (Proxy)
EXPOSE :80
EXPOSE :21
EXPOSE :8800

# Autostart script that is invoked during container start
CMD ["/usr/bin/startup"]

