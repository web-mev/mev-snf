From r-base:4.3.0

SHELL ["/bin/bash", "-c"]

# Install all software under /opt/software:
RUN mkdir -p /opt/software 

# Get the miniconda installer script and run it in silent mode:
RUN wget http://compbio.cs.toronto.edu/SNF/SNF/Software_files/SNFtool_v2.1.tar.gz -O /opt/software/SNFtool_v2.1.tar.gz && \
  cd /opt/software && \
  tar -xzf SNFtool_v2.1.tar.gz && \
  R CMD INSTALL SNFtool

ADD run_snf.R /opt/software/