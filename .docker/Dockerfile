FROM registry.fedoraproject.org/fedora-minimal:latest
MAINTAINER Jinserk Baik <jinserk.baik@gmail.com>

ENV PYTHONUNBUFFERED=1

RUN microdnf update -y \
    && microdnf install -y \
    wget bison byacc diffstat flex \
    gcc gcc-c++ gcc-gfortran \
    gettext git indent intltool libtool patch patchutils \
    redhat-rpm-config rpm-build which \
    environment-modules \
    openmpi \
    nwchem nwchem-openmpi \
    python3-devel \
    openblas-devel openblas-openmp openblas-threads \
    libxc-devel fftw-devel fftw-openmpi-devel libomp-devel \
    elpa-openmpi-devel blacs-openmpi-devel hdf5-openmpi-devel \
    && microdnf clean all
    
ARG ARCH=x86_64
ARG NPROC=2
ARG VDWXC_VER=0.4.0
ARG ELPA_VER=2021.05.002

ENV C_INCLUDE_PATH=${C_INCLUDE_PATH}:/usr/include/elpa_openmp-${ELPA_VER}
ENV LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/local/lib
ENV RUNPATH=${RUNPATH}:/usr/lib64/openmpi/lib:/usr/local/lib
   
RUN echo -e ". /etc/profile\nmodule load mpi/openmpi-${ARCH}" > /root/.bash_profile
RUN cp /root/.bash_profile /root/.zprofile
RUN . /etc/profile \
    && module load mpi/openmpi-${ARCH} \
    && pip3 install -U --no-cache-dir pip setuptools wheel pyyaml \
       pytest numpy scipy ipython jupyter python-dotenv loguru
    
WORKDIR /root
RUN . /etc/profile \
    && module load mpi/openmpi-${ARCH} \
    && wget https://launchpad.net/libvdwxc/stable/${VDWXC_VER}/+download/libvdwxc-${VDWXC_VER}.tar.gz \
    && tar -xvf libvdwxc-${VDWXC_VER}.tar.gz \
    && cd libvdwxc-${VDWXC_VER}/ \
    && ./configure CFLAGS="-fPIC" FCFLAGS="-fPIC" --with-mpi --prefix=/usr/local \
    && make -j 2 && make install \
    && cd .. \
    && rm -rf libvdwxc-${VDWXC_VER}*
    
WORKDIR /
COPY .settings/gpaw-siteconfig.py /root/.gpaw/siteconfig.py
RUN . /etc/profile \
    && module load mpi/openmpi-${ARCH} \
    && pip3 install -U --no-cache-dir matplotlib ase rdkit-pypi gpaw networkx \
    && gpaw install-data --register /root/.gpaw
    
ENV ORCA_PATH=/opt/orca501

ADD ./.downloads/orca_5_0_1_linux_x86-64_shared_openmpi411.tar.xz /opt
RUN mv /opt/orca_5_0_1_linux_x86-64_shared_openmpi411 /opt/orca501
ENV PATH="${ORCA_PATH}:${PATH}"
ENV LD_LIBRARY_PATH="${ORCA_PATH}:${LD_LIBRARY_PATH}"

ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
ENV OMP_NUM_THREADS=${NPROC}
COPY .settings/mca-params.conf /root/.openmpi/mca-params.conf

