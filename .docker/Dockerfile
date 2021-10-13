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

ADD ./.downloads/orca_5_0_1_linux_x86-64_shared_openmpi411.tar.xz /opt
RUN mv /opt/orca_5_0_1_linux_x86-64_shared_openmpi411 /opt/orca501
ENV PATH="/opt/orca501:${PATH}"
ENV LD_LIBRARY_PATH="/opt/orca501:${LD_LIBRARY_PATH}"

COPY ./.settings/openmpi-default-hostfile /etc/openmpi/
ENV ASE_NWCHEM_COMMAND="mpirun --use-hwthread-cpus nwchem PREFIX.inp >> PREFIX.out 2> PREFIX.err"
ENV ASE_ORCA_COMMAND="/opt/orca501/orca PREFIX.inp >> PREFIX.out 2> PREFIX.err"

