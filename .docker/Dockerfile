FROM python:3.9-slim AS back
ENV DEBIAN_FRONTEND noninteractive
ENV PYTHONUNBUFFERED=1
ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

RUN apt-get update -y && apt-get -y install apt-utils && apt-get -y dist-upgrade
RUN apt-get install -y --no-install-recommends \
    build-essential gfortran openbabel nwchem \
    libopenblas-dev libopenmpi-dev libxc-dev libfftw3-dev \
    libomp-dev libscalapack-openmpi-dev \
    && apt-get -y autoclean \
    && apt-get -y autoremove --purge \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install -U --no-cache-dir pip setuptools wheel pyyaml numpy scipy ipython jupyter

COPY ./.settings/gpaw-siteconfig.py /root/.gpaw/siteconfig.py
RUN pip3 install -U --no-cache-dir gpaw
RUN gpaw install-data --register /root/.gpaw
