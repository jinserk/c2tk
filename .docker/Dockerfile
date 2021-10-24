FROM registry.fedoraproject.org/fedora-minimal:latest
MAINTAINER Jinserk Baik <jinserk.baik@gmail.com>

ENV PYTHONUNBUFFERED=1

## todo with root
RUN microdnf update -y \
    && microdnf install -y \
    sudo curl wget bison byacc diffstat flex \
    zlib-devel bzip2 bzip2-devel readline-devel libffi-devel \
    sqlite sqlite-devel openssl openssl-devel \
    gcc gcc-c++ gcc-gfortran \
    gettext make git indent intltool libtool patch patchutils \
    redhat-rpm-config rpm-build which \
    environment-modules \
    openmpi \
    nwchem nwchem-openmpi \
    python3-devel python3-tkinter \
    openblas-devel openblas-openmp openblas-threads \
    libxc-devel fftw-devel fftw-openmpi-devel libomp-devel \
    elpa-openmpi-devel blacs-openmpi-devel hdf5-openmpi-devel \
    && microdnf clean all

ARG ARCH=x86_64
ARG NPROC=2
ARG USRNAME=mpiuser
ARG USRID=1000
ARG GRPID=1000
ARG VDWXC_VER=0.4.0
ARG ELPA_VER=2021.05.002

RUN groupadd -g ${GRPID} ${USRNAME} \
    && useradd -u ${USRID} -g ${GRPID} -G wheel -m -s /bin/bash ${USRNAME} \
    && sed -i /etc/sudoers -re 's/^%sudo.*/%sudo ALL=(ALL:ALL) NOPASSWD: ALL/g' \
    && sed -i /etc/sudoers -re 's/^root.*/root ALL=(ALL:ALL) NOPASSWD: ALL/g' \
    && sed -i /etc/sudoers -re 's/^#includedir.*/## **Removed the include directive** ##"/g' \
    && echo "${USRNAME} ALL=(ALL) NOPASSWD: ALL" >> /etc/sudoers \
    && echo "Customized the sudoers file for passwordless access to the ${USRNAME} user!" \
    && su - ${USRNAME} -c id

WORKDIR /root

ENV C_INCLUDE_PATH=${C_INCLUDE_PATH}:/usr/include/elpa_openmp-${ELPA_VER}
ENV LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/local/lib
ENV RUNPATH=${RUNPATH}:/usr/lib64/openmpi/lib:/usr/local/lib

RUN . /etc/profile \
    && module load mpi/openmpi-${ARCH} \
    && wget https://launchpad.net/libvdwxc/stable/${VDWXC_VER}/+download/libvdwxc-${VDWXC_VER}.tar.gz \
    && tar -xvf libvdwxc-${VDWXC_VER}.tar.gz \
    && cd libvdwxc-${VDWXC_VER}/ \
    && ./configure CFLAGS="-fPIC" FCFLAGS="-fPIC" --with-mpi --prefix=/usr/local \
    && make -j \
    && make install \
    && cd .. \
    && rm -rf libvdwxc-${VDWXC_VER}*

RUN mkdir -p /root/.openmpi
COPY .settings/mca-params.conf /root/.openmpi/mca-params.conf
ENV OMP_NUM_THREADS=${NPROC}

COPY .settings/.bash_profile /root
ENV C_INCLUDE_PATH=${C_INCLUDE_PATH}:/usr/include/elpa_openmp-${ELPA_VER}
ENV LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/local/lib

RUN . /root/.bash_profile \
    && pip install -U --no-cache-dir pip setuptools wheel pyyaml \
       pytest numpy scipy mpi4py ipython jupyter python-dotenv loguru

RUN mkdir -p /root/.gpaw
COPY .settings/gpaw-siteconfig.py /root/.gpaw/siteconfig.py
RUN . /root/.bash_profile \
    && pip install -U --no-cache-dir matplotlib networkx \
       ase rdkit-pypi gpaw nglview plams matscipy pyamg

ENV ORCA_PATH=/opt/orca501
ADD .downloads/orca_5_0_1_linux_x86-64_shared_openmpi411.tar.xz /opt
RUN mv /opt/orca_5_0_1_linux_x86-64_shared_openmpi411 ${ORCA_PATH}

## todo with user
USER ${USRNAME}
ENV HOME=/home/${USRNAME}
WORKDIR ${HOME}

RUN mkdir -p ${HOME}/.openmpi
COPY .settings/mca-params.conf ${HOME}/.openmpi/mca-params.conf
ENV OMP_NUM_THREADS=${NPROC}

COPY .settings/.bash_profile ${HOME}
ENV C_INCLUDE_PATH=${C_INCLUDE_PATH}:/usr/include/elpa_openmp-${ELPA_VER}
ENV LD_LIBRARY_PATH="${ORCA_PATH}:${LD_LIBRARY_PATH}"
ENV RUNPATH=${RUNPATH}:/usr/lib64/openmpi/lib:/usr/local/lib
ENV PATH=${HOME}/.local/bin:${PATH}

ENV ORCA_PATH=/opt/orca501
ENV PATH="${ORCA_PATH}:${PATH}"

RUN mkdir -p ${HOME}/.ipython/profile_default
COPY .settings/ipython_config.py ${HOME}/.ipython/profile_default/

RUN mkdir -p ${HOME}/.jupyter
COPY .settings/jupyter_notebook_config.py ${HOME}/.jupyter/

RUN mkdir -p ${HOME}/.gpaw
RUN . ${HOME}/.bash_profile \
    && gpaw install-data --register ${HOME}/.gpaw \
    && jupyter-nbextension enable nglview --py --user
