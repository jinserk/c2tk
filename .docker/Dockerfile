FROM python:3.9-slim AS back

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update -y && apt-get -y install apt-utils && apt-get -y dist-upgrade
RUN apt-get install -y --no-install-recommends \
    wget sudo tcsh build-essential gawk libtool gcc gfortran \
    mpi-default-bin mpi-default-dev libscalapack-mpi-dev libomp-dev \
    libopenblas-dev libopenblas64-dev libopenblas-openmp-dev libopenblas64-openmp-dev \
    libxc-dev libfftw3-dev libfftw3-mpi-dev libhdf5-dev \
    openbabel psi4 psi4-data \
    && apt-get -y autoclean \
    && apt-get -y autoremove --purge \
    && rm -rf /var/lib/apt/lists/*

ARG NPROC=2
ARG USERNAME=mpiuser
ARG USERID=1000
ARG GROUPID=1000

ARG VDWXC_VER=0.4.0
ARG ELPA_VER=2021.11.001

RUN groupadd -f -g ${GROUPID} ${USERNAME} \
  && useradd -u ${USERID} -g ${GROUPID} -G sudo -m -s /bin/bash ${USERNAME} \
  && echo "${USERNAME} ALL=(ALL:ALL) NOPASSWD:ALL" > /etc/sudoers.d/${USERNAME}
  
USER ${USERNAME}
ENV HOME=/home/${USERNAME}
WORKDIR ${HOME}

RUN wget https://launchpad.net/libvdwxc/stable/${VDWXC_VER}/+download/libvdwxc-${VDWXC_VER}.tar.gz \
    && tar -xvf libvdwxc-${VDWXC_VER}.tar.gz \
    && cd libvdwxc-${VDWXC_VER}/ \
    && ./configure CFLAGS="-fPIC" FCFLAGS="-fPIC" --with-mpi --prefix=/usr/local \
    && make -j \
    && sudo make install \
    && cd .. \
    && rm -rf libvdwxc-${VDWXC_VER}*

RUN wget https://elpa.mpcdf.mpg.de/software/tarball-archive/Releases/${ELPA_VER}/elpa-${ELPA_VER}.tar.gz \
    && tar -xvf elpa-${ELPA_VER}.tar.gz \
    && cd elpa-${ELPA_VER}/ \
    && ./configure CC="mpicc" FC="mpif90" AWK="gawk" \
        FCFLAGS="-O2 -mmmx -msse -msse2 -msse3 -mssse3 -msse4.1 -msse4.2 -maes -mavx -mfma -mavx2" \
        CFLAGS="-O2 -mmmx -msse -msse2 -msse3 -mssse3 -msse4.1 -msse4.2 -maes -mavx -mfma -mavx2" \
        --enable-openmp --with-mpi=yes --disable-avx512 --prefix=/usr/local \
    && make -j \
    && sudo make install \
    && cd .. \
    && rm -rf elpa-${ELPA_VER}*
    
ENV PATH=${HOME}/.local/bin:${PATH}
ENV PYTHONUNBUFFERED=1
ENV OMP_NUM_THREADS=${NPROC}
#ENV OMPI_ALLOW_RUN_AS_ROOT=1
#ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
#ENV C_INCLUDE_PATH=${C_INCLUDE_PATH}:/usr/local/include/elpa-${ELPA_VER}
ENV LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/local/lib
ENV PYTHONPATH=${PYTHONPATH}:/usr/lib/x86_64-linux-gnu

RUN mkdir -p ${HOME}/.pip
COPY --chown=${USERID}:${GROUPID} .settings/pip.conf ${HOME}/.pip/
RUN mkdir -p ${HOME}/.openmpi
COPY --chown=${USERID}:${GROUPID} .settings/mca-params.conf ${HOME}/.openmpi/mca-params.conf

RUN pip install -U pip setuptools wheel pyyaml pytest
RUN pip install -U pybind11 mpi4py cython ipython python-dotenv loguru
RUN pip install -U numpy scipy h5py jupyterlab matplotlib networkx
 
RUN mkdir -p ${HOME}/.gpaw
COPY --chown=${USERID}:${GROUPID} .settings/gpaw-siteconfig.py ${HOME}/.gpaw/siteconfig.py
COPY --chown=${USERID}:${GROUPID} ./requirements.txt /tmp/requirements.txt
RUN pip install -U -r /tmp/requirements.txt
RUN gpaw install-data --register ${HOME}/.gpaw

RUN mkdir -p ${HOME}/.ipython/profile_default
COPY --chown=${USERID}:${GROUPID} .settings/ipython_config.py ${HOME}/.ipython/profile_default/

RUN mkdir -p ${HOME}/.jupyter/lab/user-settings
#COPY .settings/jupyter_notebook_config.py ${HOME}/.jupyter/
COPY --chown=${USERID}:${GROUPID} .settings/jupyter_server_config.py ${HOME}/.jupyter/

ENV ORCA_PATH=/opt/orca502
ENV PATH="${ORCA_PATH}:${PATH}"
ENV LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${ORCA_PATH}
ADD --chown=root:root ./.downloads/orca_5_0_2_linux_x86-64_shared_openmpi411.tar.xz /opt
RUN sudo mv /opt/orca_5_0_2_linux_x86-64_shared_openmpi411 ${ORCA_PATH}

#ENV GAMESS_PATH=/opt/gamess
#ADD --chown=${USERID}:${GROUPID} ./.downloads/gamess-current.tar.gz /opt
#RUN cd ${GAMESS_PATH} \
#    && configure --help
