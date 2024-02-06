FROM mambaorg/micromamba

# for building the conda env
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
# Build the environment based on the yml file
RUN micromamba install -y -n base -f /tmp/environment.yml \
  && micromamba clean --all --yes

USER $MAMBA_USER
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN curl -o /tmp/snf.tar.gz https://cran.r-project.org/src/contrib/SNFtool_2.3.1.tar.gz && \
  cd /tmp && \
  tar -xzf snf.tar.gz && \
  R CMD INSTALL SNFtool

ADD run_snf.R /usr/local/bin/

