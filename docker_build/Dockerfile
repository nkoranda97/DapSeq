FROM condaforge/mambaforge:latest

RUN conda config --set remote_max_retries 10 \
 && conda config --set remote_backoff_factor 5 \
 && conda config --set remote_connect_timeout_secs 60 \
 && conda config --set remote_read_timeout_secs 120

COPY environment.yml .
RUN mamba env create -f environment.yml && mamba clean -afy
ENV PATH=/opt/conda/envs/dapseq/bin:$PATH
