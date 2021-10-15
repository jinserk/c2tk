FROM registry.fedoraproject.org/fedora-minimal:latest
ENV PYTHONUNBUFFERED=1

RUN microdnf update -y \
    && microdnf install -y \
    wget bison byacc diffstat flex \
    gcc gcc-c++ gcc-gfortran \
    gettext git indent intltool libtool patch patchutils \
    redhat-rpm-config rpm-build which \
    environment-modules \
    mpich \
    nwchem nwchem-mpich \
    python3-devel \
    openblas-devel openblas-openmp openblas-threads \
    libxc-devel fftw-devel fftw-mpich-devel libomp-devel \
    elpa-mpich-devel blacs-mpich-devel hdf5-mpich-devel \
    && microdnf clean all
    
ARG ARCH=x86_64
ARG VDWXC_VER=0.4.0
ARG ELPA_VER=2021.05.002

ENV C_INCLUDE_PATH=${C_INCLUDE_PATH}:/usr/include/elpa_openmp-${ELPA_VER}
ENV LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/local/lib
ENV RUNPATH=${RUNPATH}:/usr/lib64/mpich/lib:/usr/local/lib
   
RUN echo -e ". /etc/profile\nmodule load mpi/mpich-${ARCH}" > /root/.bash_profile
RUN cp /root/.bash_profile /root/.zprofile
RUN . /etc/profile \
    && module load mpi/mpich-${ARCH} \
    && pip3 install -U --no-cache-dir pip setuptools wheel pyyaml pytest numpy scipy ipython jupyter
    
WORKDIR /root
RUN . /etc/profile \
    && module load mpi/mpich-${ARCH} \
    && wget https://launchpad.net/libvdwxc/stable/${VDWXC_VER}/+download/libvdwxc-${VDWXC_VER}.tar.gz \
    && tar -xvf libvdwxc-${VDWXC_VER}.tar.gz \
    && cd libvdwxc-${VDWXC_VER}/ \
    && ./configure CFLAGS="-fPIC" FCFLAGS="-fPIC" --with-mpi --prefix=/usr/local \
    && make -j 2 && make install \
    && cd .. \
    && rm -rf libvdwxc-${VDWXC_VER}*
    
WORKDIR /
COPY ./.settings/gpaw-siteconfig.py /root/.gpaw/siteconfig.py
RUN . /etc/profile \
    && module load mpi/mpich-${ARCH} \
    && pip3 install -U --no-cache-dir matplotlib ase rdkit-pypi gpaw networkx \
    && gpaw install-data --register /root/.gpaw
    
#ADD ./.downloads/orca_5_0_1_linux_x86-64_shared_openmpi411.tar.xz /opt
#RUN mv /opt/orca_5_0_1_linux_x86-64_shared_openmpi411 /opt/orca501
#ENV PATH="/opt/orca501:${PATH}"
#ENV LD_LIBRARY_PATH="/opt/orca501:${LD_LIBRARY_PATH}"

ARG NPROC=1
ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

ENV ASE_NWCHEM_COMMAND="mpirun -np ${NPROC} nwchem_mpich PREFIX.nwi > PREFIX.nwo"
#ENV ASE_ORCA_COMMAND="/opt/orca501/orca PREFIX.inp >> PREFIX.out 2> PREFIX.err"

