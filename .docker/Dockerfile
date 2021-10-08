FROM python:3.9-slim AS back
ENV DEBIAN_FRONTEND noninteractive
ENV PYTHONUNBUFFERED=1

RUN apt-get update -y && apt-get -y install apt-utils && apt-get -y dist-upgrade
RUN apt-get install -y --no-install-recommends \
    openbabel nwchem \
    && apt-get -y autoclean \
    && apt-get -y autoremove --purge \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install -U --no-cache-dir pip setuptools wheel pyyaml numpy scipy ipython jupyter

