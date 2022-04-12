############
# Build base
from ubuntu:18.04 as build-base
RUN apt-get update
RUN apt-get install -y \
      build-essential \
      git \
      wget \
      zlib1g-dev


###########
# Build BWA
from build-base as build-bwa
RUN apt-get install -y \
      build-essential \
      git \
      g++ \
      gcc \
      git \
      openjdk-8-jre \
      perl \
      python \
      zlib1g-dev
# Install BWA
RUN   git clone https://github.com/lh3/bwa.git \
      && cd bwa \
      && git checkout 0.7.12 \
      && make

################
# Build samtools
FROM build-base as build-samtools
RUN apt-get install -y \
    build-essential \
    autoconf \
    g++ \
    git \
    libbz2-dev \
    liblzma-dev \
    make \
    ncurses-dev \
    wget \
    zlib1g-dev
RUN wget --quiet https://github.com/samtools/samtools/releases/download/1.7/samtools-1.7.tar.bz2 \
      && tar xf samtools-1.7.tar.bz2 \
      && cd samtools-1.7 \
      && make install

#####################
# Build trim-adapters
FROM build-base as build-trim-adapters
RUN apt-get install -y \
      build-essential \
      libboost-dev \
      git \
      zlib1g-dev
RUN git clone https://bitbucket.org/jvierstra/bio-tools.git \
      && cd bio-tools \
      && git checkout 6fe54fa5a3 \
      && make

########
# Picard
from build-base as get-picard
RUN wget --quiet https://github.com/broadinstitute/picard/releases/download/2.8.1/picard.jar

#######
# fastp
from build-base as get-fastp
RUN wget --quiet http://opengene.org/fastp/fastp.0.23.2 \
    && mv fastp.0.23.2 fastp \
    && chmod ugo+x fastp

########
# Bedops
from build-base as build-bedops
RUN apt-get install -y \
      build-essential \
      git \
      libbz2-dev
RUN git clone https://github.com/bedops/bedops.git \
      && cd bedops \
      && git checkout v2.4.35 \
      && make \
      && make install

##########
# Hotspot1
from build-base as build-hotspot1
RUN apt-get install -y \
      build-essential \
      git \
      libgsl-dev \
      wget
RUN git clone https://github.com/StamLab/hotspot.git \
      && cd hotspot \
      && git checkout v4.1.1 \
      && cd hotspot-distr/hotspot-deploy \
      && make

###########
# Kentutils
from build-base as build-kentutils
RUN apt-get install -y \
      build-essential \
      git \
      libmysqlclient-dev \
      libpng-dev \
      libssh-dev \
      wget \
      zlib1g-dev
RUN wget --quiet https://github.com/ENCODE-DCC/kentUtils/archive/v302.0.0.tar.gz \
      && tar xf v302.0.0.tar.gz \
      && cd kentUtils-302.0.0 \
      && make

##########
# Bedtools
from build-base as build-bedtools
RUN apt-get install -y \
      python
RUN wget --quiet https://github.com/arq5x/bedtools2/releases/download/v2.25.0/bedtools-2.25.0.tar.gz \
      && tar xf bedtools-2.25.0.tar.gz \
      && cd bedtools2 \
      && make

########
# Preseq
from build-base as build-preseq
RUN apt-get install -y \
      libgsl-dev
RUN git clone --recurse-submodules https://github.com/smithlabcode/preseq.git \
   && cd preseq \
   && git checkout v2.0.1 \
   && git submodule update \
   && make

##########
# Hotspot2
from build-base as build-hotspot2
RUN git clone https://github.com/Altius/hotspot2.git \
  && cd hotspot2 \
  && make \
  && cd / \
  && git clone https://github.com/StamLab/modwt.git \
  && cd modwt \
  && git checkout 28e9f479c737836ffc870199f2468e30659ab38d \
  && make


#######################
# Final image for DNase
from ubuntu:18.04 as stampipes-dnase

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update
RUN apt-get install -y \
      bash \
      bc \
      bowtie \
      build-essential \
      coreutils \
      gawk \
      libboost-dev \
      libgsl-dev \
      littler \
      openjdk-8-jre \
      python-dev \
      python-pip \
      python3 \
      python3-pip \
      tabix \
      wget \
      zlib1g-dev

COPY ./requirements.pip.txt /stampipes/
RUN pip install -r /stampipes/requirements.pip.txt
RUN pip3 install -r /stampipes/requirements.pip.txt

COPY ./scripts /stampipes/scripts
COPY ./processes /stampipes/processes
COPY ./makefiles /stampipes/makefiles
COPY ./awk /stampipes/awk
ENV STAMPIPES=/stampipes

# Copy in dependencies
COPY --from=build-bwa /bwa/bwa /usr/local/bin/
COPY --from=build-trim-adapters /bio-tools/apps/trim-adapters-illumina/trim-adapters-illumina /usr/local/bin/
COPY --from=build-samtools /usr/local/bin/samtools /usr/local/bin
COPY --from=build-bedops /bedops/bin /usr/local/bin
ENV HOTSPOT_DIR /hotspot
COPY --from=build-hotspot1 /hotspot/hotspot-distr/ $HOTSPOT_DIR
COPY --from=build-hotspot2 /hotspot2/bin /usr/local/bin/
COPY --from=build-hotspot2 /hotspot2/scripts /usr/local/bin/
COPY --from=build-hotspot2 /modwt/bin /usr/local/bin/
COPY --from=build-kentutils /kentUtils-302.0.0/bin/ /usr/local/bin/
COPY --from=build-bedtools /bedtools2/bin/ /usr/local/bin/
COPY --from=build-preseq /preseq/preseq /usr/local/bin/
COPY --from=get-picard /picard.jar /usr/local/lib/picard.jar
COPY --from=get-fastp /fastp /usr/local/bin/

# Make alias for picard
RUN echo -e '#!/bin/bash\njava -jar /usr/local/lib/picard.jar $@' \
      > /usr/local/bin/picard \
      && chmod +x /usr/local/bin/picard

##########
# RNA deps
##########

######
# STAR
from build-base as build-star
ARG star_version=2.4.2a
RUN wget --quiet https://github.com/alexdobin/STAR/archive/refs/tags/STAR_${star_version}.tar.gz \
    && tar -xf STAR_${star_version}.tar.gz \
    && cp /STAR-STAR_${star_version}/bin/Linux_x86_64/* /usr/local/bin/

######
# RSEM
FROM build-base as build-rsem
ARG rsem_version=1.2.30
RUN wget --quiet https://github.com/deweylab/RSEM/archive/refs/tags/v${rsem_version}.tar.gz \
    && tar -xf v${rsem_version}.tar.gz \
    && find RSEM-${rsem_version} -maxdepth 1 -executable -type f -exec cp '{}' /usr/local/bin ';'

###########
# Cufflinks
FROM build-base as build-cufflinks
ARG cufflinks_version=2.2.1
RUN wget --quiet http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-${cufflinks_version}.Linux_x86_64.tar.gz \
    && tar -xf cufflinks-${cufflinks_version}.Linux_x64_64.tar.gz \
    && find cufflinks-${cufflinks_version}.Linux_x86_64 -maxdepth 1 -executable -type f -exec cp '{}' /usr/local/bin ';'

#########
# Anaquin
FROM build-base as build-anaquin
ARG anaquin_version=2.0.1
# TODO - where do we find this??

###########
# Stringtie
FROM build-base as build-stringtie
ARG stringtie_version=1.3.4d
RUN wget --quiet http://ccb.jhu.edu/software/stringtie/dl/stringtie-${stringtie_version}.Linux_x86_64.tar.gz \
    && tar -xf stringtie-${stringtie_version}.Linux_x86_64.tar.gz \
    && cp stringtie-${stringtie_version}.Linux_x86_64/stringtie /usr/local/bin/

#########
# Subread
FROM build-base as build-subread
ARG subread_version=1.5.1
RUN wget --quiet https://downloads.sourceforge.net/project/subread/subread-${subread_version}/subread-${subread_version}-Linux-x86_64.tar.gz \
    && tar -xf subread-${subread_version}-Linux-x86_64.tar.gz \
    && find subread-${subread_version}-Linux-x86_64/bin -executable -type f -exec cp '{}' /usr/local/bin ';'

##########
# Kallisto
FROM build-base as build-kallisto
ARG kallisto_version=0.43.1
ARG hdf5_version=1.10.1
RUN apt-get -y install cmake
RUN wget --quiet https://support.hdfgroup.org/ftp/HDF5/prev-releases/hdf5-1.10/hdf5-${hdf5_version}/src/hdf5-${hdf5_version}.tar.bz2 \
      && tar xf hdf5-${hdf5_version}.tar.bz2 \
      && cd  hdf5-${hdf5_version} \
      && ./configure --prefix / \
      && make install
RUN wget --quiet https://github.com/pachterlab/kallisto/archive/v${kallisto_version}.tar.gz \
      && tar xf v${kallisto_version}.tar.gz \
      && cd kallisto-${kallisto_version} \
      && mkdir build \
      && cd build \
      && cmake .. \
      && make install

#####################
# Final RNA-seq image
from ubuntu:18.04 as stampipes-rna
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update
RUN apt-get install -y \
      bash \
      bc \
      bowtie \
      build-essential \
      coreutils \
      gawk \
      libboost-dev \
      libgsl-dev \
      littler \
      openjdk-8-jre \
      python-dev \
      python-pip \
      python3 \
      python3-pip \
      tabix \
      wget \
      zlib1g-dev

COPY ./requirements.pip.txt /stampipes/
RUN pip install -r /stampipes/requirements.pip.txt
RUN pip3 install -r /stampipes/requirements.pip.txt

COPY ./scripts /stampipes/scripts
COPY ./processes /stampipes/processes
COPY ./makefiles /stampipes/makefiles
COPY ./awk /stampipes/awk
ENV STAMPIPES=/stampipes

# Copy in dependencies
COPY --from=build-trim-adapters /bio-tools/apps/trim-adapters-illumina/trim-adapters-illumina /usr/local/bin/
COPY --from=build-samtools /usr/local/bin/samtools /usr/local/bin
COPY --from=build-bedops /bedops/bin /usr/local/bin
COPY --from=build-kentutils /kentUtils-302.0.0/bin/ /usr/local/bin/
COPY --from=build-bedtools /bedtools2/bin/ /usr/local/bin/
COPY --from=build-preseq /preseq/preseq /usr/local/bin/
COPY --from=get-picard /picard.jar /usr/local/lib/picard.jar
COPY --from=get-fastp /fastp /usr/local/bin/
COPY --from=build-star /usr/local/bin/STAR* /usr/local/bin/
COPY --from=build-rsem /usr/local/bin/* /usr/local/bin/
COPY --from=build-kallisto /usr/local/bin/kallisto /usr/local/bin/
COPY --from=build-kallisto /lib/libhdf5.so.101 /lib
COPY --from=build-subread /usr/local/bin/* /usr/local/bin/
COPY --from=build-stringtie /usr/local/bin/stringtie /usr/local/bin/

# Make alias for picard
RUN echo -e '#!/bin/bash\njava -jar /usr/local/lib/picard.jar $@' \
      > /usr/local/bin/picard \
      && chmod +x /usr/local/bin/picard
