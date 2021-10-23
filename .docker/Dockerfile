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
    
## todo with user
USER ${USRNAME}
ENV HOME=/home/${USRNAME}
WORKDIR ${HOME}

RUN mkdir -p ${HOME}/.openmpi
COPY .settings/mca-params.conf ${HOME}/.openmpi/mca-params.conf
ENV OMP_NUM_THREADS=${NPROC}
 
ENV C_INCLUDE_PATH=${C_INCLUDE_PATH}:/usr/include/elpa_openmp-${ELPA_VER}
ENV LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/local/lib
ENV RUNPATH=${RUNPATH}:/usr/lib64/openmpi/lib:/usr/local/lib
ENV PATH=${HOME}/.local/bin:${PATH}
    
COPY .settings/.bash_profile ${HOME}
COPY .settings/.bash_profile ${HOME}/.profile
COPY .settings/.bashrc ${HOME}

RUN git clone git://github.com/pyenv/pyenv.git .pyenv 

ENV PYENV_ROOT=${HOME}/.pyenv
ENV PATH=${PYENV_ROOT}/shims:${PYENV_ROOT}/bin:${PATH}
        
RUN pyenv install 3.9.7 \
    && pyenv global 3.9.7 \
    && pyenv rehash

RUN . ${HOME}/.bash_profile \
    && pip install -U --no-cache-dir pip setuptools wheel pyyaml \
       pytest numpy scipy mpi4py ipython jupyter python-dotenv loguru
   
RUN mkdir -p ${HOME}/.gpaw
COPY .settings/gpaw-siteconfig.py ${HOME}/.gpaw/siteconfig.py

RUN . ${HOME}/.bash_profile \
    && pip install -U --no-cache-dir matplotlib networkx \
       ase rdkit-pypi gpaw nglview plams matscipy pyamg \
    && gpaw install-data --register ${HOME}/.gpaw
    
RUN jupyter-nbextension enable nglview --py --user

ENV ORCA_PATH=/opt/orca501

ADD ./.downloads/orca_5_0_1_linux_x86-64_shared_openmpi411.tar.xz ${HOME}
RUN sudo mv ${HOME}/orca_5_0_1_linux_x86-64_shared_openmpi411 ${ORCA_PATH}

ENV PATH="${ORCA_PATH}:${PATH}"
ENV LD_LIBRARY_PATH="${ORCA_PATH}:${LD_LIBRARY_PATH}"
