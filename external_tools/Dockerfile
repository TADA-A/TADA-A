FROM debian:testing-slim

WORKDIR /data
ENV MINICONDA_VERSION 4.4.10
ENV PATH /opt/miniconda3/bin:$PATH

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    curl ca-certificates bzip2 unzip \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log

# Install Python packages
RUN curl https://repo.continuum.io/miniconda/Miniconda3-$MINICONDA_VERSION-Linux-x86_64.sh -o MCON.sh \
    && /bin/bash MCON.sh -b -p /opt/miniconda3 \
    && ln -s /opt/miniconda3/etc/profile.d/conda.sh /etc/profile.d/conda.sh \
    && conda install -c conda-forge python-xxhash psutil \
    && conda clean --all -tipsy && rm -rf /tmp/* $HOME/.cache

RUN pip install sos && rm -rf $HOME/.cache

RUN curl -L http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa \
    -o /usr/local/bin/twoBitToFa

RUN curl -L http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig \
    -o /usr/local/bin/wigToBigWig

COPY base_mutarate_to_wiggle_file.sh /usr/local/bin/
COPY MutRateBase_for_TADA-A_v2.py /usr/local/bin/
COPY tri_extract_for_TADA-A.py /usr/local/bin/
COPY bedtools-2.17.0/bin/* /usr/local/bin/
COPY fordist_1KG_mutation_rate_table.txt ./
COPY hg19.genome.txt ./
RUN chmod 755 /usr/local/bin/*

cmd ["bash"]