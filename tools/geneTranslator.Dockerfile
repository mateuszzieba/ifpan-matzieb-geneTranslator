# geneTranslator.Dockerfile

FROM rocker/tidyverse:latest

WORKDIR /home/rstudio/project

COPY . .

# Zainstaluj BiocManager
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org/')"

# Zainstaluj pakiety Bioconductor i CRAN
RUN R -e "BiocManager::install(c(\
  'biomaRt', \
  'HGNChelper', \
  'UniProt.ws', \
  'org.Hs.eg.db' \
), ask=FALSE, update=FALSE)" \
 && install2.r --error rentrez

CMD ["R"]
