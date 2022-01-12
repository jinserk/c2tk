FROM archlinux/archlinux:base-devel
MAINTAINER Jinserk Baik <jinserk.baik@gmail.com>

RUN pacman-key --init
RUN pacman -Syyu --needed --noconfirm \
    && pacman -Syu --needed --noconfirm base-devel git go gcc-fortran \
    && pacman -Scc --noconfirm

# create user and group
ARG USERNAME=mpiuser
ARG USERID=1000
ARG GROUPID=1000
ARG NPROC=1

RUN groupadd -f -g ${GROUPID} ${USERNAME} \
  && useradd -u ${USERID} -g ${GROUPID} -G wheel -m -s /bin/bash ${USERNAME} \
  && echo "${USERNAME} ALL=(ALL:ALL) NOPASSWD:ALL" > /etc/sudoers.d/${USERNAME}
  
USER ${USERNAME}
ENV HOME=/home/${USERNAME}
WORKDIR ${HOME}

RUN git clone https://aur.archlinux.org/yay.git \
    && cd yay \
    && makepkg -sri --needed --noconfirm \
    && cd .. && rm -rf yay \
    && yay -Syu --devel --needed --noconfirm \
    && yay -Scc --noconfirm

RUN yay -Syu --needed --noconfirm \
    python-pip openmpi openmp openblas-lapack fftw elpa libvdwxc libxc \
    && yay -Scc --noconfirm

ENV PYTHONUNBUFFERED=1
ENV PATH=${HOME}/.local/bin:${PATH}
#ENV PYTHONPATH=${HOME}/.local/lib/python3.10/site-packages:/usr/lib/python3/dist-packages:${PYTHONPATH}
ENV OMP_NUM_THREADS=${NPROC}

RUN mkdir -p ${HOME}/.pip
COPY --chown=${USERID}:${GROUPID} .settings/pip.conf ${HOME}/.pip/

RUN pip install -U pip setuptools wheel pyyaml pytest
RUN pip install -U pybind11 ipython python-dotenv loguru
RUN pip install -U numpy scipy mpi4py jupyterlab matplotlib networkx
    
RUN mkdir -p ${HOME}/.gpaw
COPY --chown=${USERID}:${GROUPID} .settings/gpaw-siteconfig.py ${HOME}/.gpaw/siteconfig.py
COPY --chown=${USERID}:${GROUPID} ./requirements.txt /tmp/requirements.txt
RUN pip install -U -r /tmp/requirements.txt

RUN mkdir -p ${HOME}/.openmpi
COPY --chown=${USERID}:${GROUPID} .settings/mca-params.conf ${HOME}/.openmpi/mca-params.conf
#RUN gpaw install-data --register ${HOME}/.gpaw \
#    && jupyter-nbextension enable nglview --py --user
    
RUN mkdir -p ${HOME}/.ipython/profile_default
COPY --chown=${USERID}:${GROUPID} .settings/ipython_config.py ${HOME}/.ipython/profile_default/

RUN mkdir -p ${HOME}/.jupyter/lab/user-settings
#COPY .settings/jupyter_notebook_config.py ${HOME}/.jupyter/
COPY --chown=${USERID}:${GROUPID} .settings/jupyter_server_config.py ${HOME}/.jupyter/

