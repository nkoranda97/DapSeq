FROM condaforge/mambaforge:latest

RUN conda config --set remote_max_retries 10 \
 && conda config --set remote_backoff_factor 5 \
 && conda config --set remote_connect_timeout_secs 60 \
 && conda config --set remote_read_timeout_secs 120

COPY environment.yml .
RUN mamba env create -f environment.yml && mamba clean -afy
ENV PATH=/opt/conda/envs/dapseq/bin:$PATH

# Install GEM peak caller
RUN mkdir -p /opt/gem \
 && wget -q http://groups.csail.mit.edu/cgs/gem/download/gem.v3.4.tar.gz \
    -O /tmp/gem.tar.gz \
 && tar -xzf /tmp/gem.tar.gz -C /opt/gem --strip-components=1 \
 && rm /tmp/gem.tar.gz