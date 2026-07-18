FROM docker.io/mambaorg/micromamba:2.6.2

ARG MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH=/opt/conda/bin:$PATH

COPY --chown=$MAMBA_USER:$MAMBA_USER docker/environment.yml /tmp/mutt-environment.yml
RUN micromamba install --yes --name base --file /tmp/mutt-environment.yml && \
    micromamba clean --all --yes

WORKDIR /tmp/mutt
COPY --chown=$MAMBA_USER:$MAMBA_USER . /tmp/mutt
RUN bash -n /tmp/mutt/scripts/hpc_run_functional_all.sh && \
    bash -n /tmp/mutt/zip-push-gitlfs.sh && \
    Rscript --vanilla -e 'files <- list.files("/tmp/mutt/scripts", pattern = "[.]R$", full.names = TRUE); invisible(lapply(files, parse))' && \
    Rscript --vanilla -e 'testthat::test_local("/tmp/mutt", stop_on_failure = TRUE)' && \
    R CMD INSTALL --no-multiarch --with-keep.source /tmp/mutt && \
    Rscript --vanilla -e 'stopifnot(identical(getNamespaceExports("mutt"), "mutt"), dir.exists(system.file("extdata", "FAPROTAX_1.2.12", package = "mutt")))' && \
    picrust2_pipeline.py --version

WORKDIR /work
CMD ["R", "--no-save"]
