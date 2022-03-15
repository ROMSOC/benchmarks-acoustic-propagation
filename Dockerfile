# Copyright (C) 2015-2022 by the RBniCS authors
#
# This file is part of RBniCS.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

FROM ubuntu:20.04
ENV LAST_UPDATED=2022-03-01

USER root
WORKDIR /tmp

# Install FEniCS
USER root

# Install FEniCS from system packages
RUN apt-get update --yes
RUN apt-get install --yes --no-install-recommends software-properties-common
RUN apt-add-repository ppa:fenics-packages/fenics
RUN apt-get update --yes && \
  apt-get install --yes --no-install-recommends python3-pip fenics

ARG NB_USER="jovyan"
ARG NB_UID="1000"
ARG NB_GID="100"

# Fix DL4006
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# Install JupyterLab
USER root

RUN chmod 777 /tmp        && \
    DEBIAN_FRONTEND="noninteractive" && \
    pip3 install --upgrade --no-cache-dir \
      jupyterlab                    \
      jupytext                      \
      jupyter-book                  \
      ghp-import

# allow jupyterlab for ipyvtk
ENV JUPYTER_ENABLE_LAB=yes
ENV PYVISTA_USE_IPYVTK=true

# Install basic Python packages
USER root
RUN apt-get update --yes
RUN apt-get install --yes --no-install-recommends python3-dev && \
    pip3 -q install sympy numpy scipy matplotlib 

# Install Pyvista
USER root
RUN apt-get update --yes && \
    apt-get install  -yq --no-install-recommends \
      libgl1 \
      libgl1-mesa-glx xvfb \
      libfontconfig1 \
      libxrender1 \
      libosmesa6  && \
    pip3 -q install cmocean colorcet imageio-ffmpeg imageio ipygany \ 
      ipyvtklink meshio panel trimesh pythreejs pyvista \ 
      piglet pyvirtualdisplay

# Install binaries from Salome
USER root
RUN apt-get update --yes && \ 
  apt-get install -yq --no-install-recommends python3-babel python3-pytest-cython python3-jinja2 \
    python3-pil python3-pyqt5 python3-pygments python3-sphinx \ 
    python3-alabaster python3-click libcminpack1 python3-cycler \
    python3-dateutil python3-docutils fftw-dev libfreeimage3 graphviz \
    python3-imagesize python3-kiwisolver clang python3-markupsafe \
    python3-matplotlib python3-packaging python3-pandas python3-psutil \
    python3-pyparsing python3-tz libqwt-qt5-6 libdc1394-22 libegl1 libexif12 \
    libglu1-mesa libgphoto2-6 libjbig0 libncurses5 libraw1394-11 libtiff5 \
    libusb-1.0-0 libxcb-xkb1 libxft2 libxi6 libxkbcommon0 libxkbcommon-x11-0 \
    libxss1 libopenexr24 libwebp6 libilmbase24 libgphoto2-port12 libopengl0 \
    python3-scipy python3-sip python3-stemmer python3-sphinx-rtd-theme \
    python3-sphinxcontrib.websupport sphinx-intl python3-statsmodels libtbb2 libtcl8.6 libtk8.6

# Install Salome and ensure that Python3 is used in any cricunstances
USER root
RUN ln -s /usr/bin/python3 /usr/bin/python
WORKDIR /usr/local
RUN apt-get update --yes && \ 
  apt-get install -yq --no-install-recommends wget
RUN wget https://files.salome-platform.org/Salome/Salome9.8.0/SALOME-9.8.0-native-UB20.04-SRC.tar.gz
RUN tar xvzf SALOME-9.8.0-native-UB20.04-SRC.tar.gz && \
    rm SALOME-9.8.0-native-UB20.04-SRC.tar.gz
RUN ln -s /usr/local/SALOME-9.8.0-native-UB20.04-SRC/BINARIES-UB20.04/petsc \
        /usr/local/SALOME-9.8.0-native-UB20.04-SRC/BINARIES-UB20.04/petsc/arch-linux-c-opt && \
    ln -s /usr/local/SALOME-9.8.0-native-UB20.04-SRC/salome /usr/local/bin/salome

# remove unnecessary files
RUN rm -rf /usr/local/SALOME-9.8.0-native-UB20.04-SRC/ARCHIVES/ && \
    rm -rf /usr/local/SALOME-9.8.0-native-UB20.04-SRC/SOURCES/ && \
    # remove unnecessary binaries - be careful here, some binaries are required for testing!
    rm -rf $HOME/salome_dir/BINARIES*/MEDCOUPLING/ && \
    rm -rf $HOME/salome_dir/BINARIES*/MeshGems/ && \
    rm -rf $HOME/salome_dir/BINARIES*/FIELDS/ && \
    rm -rf $HOME/salome_dir/BINARIES*/YACS/ && \
    rm -rf $HOME/salome_dir/BINARIES*/HOMARD/ && \
    rm -rf $HOME/salome_dir/BINARIES*/PARAVIS/ && \
    rm -rf $HOME/salome_dir/BINARIES*/HEXABLOCK/ 

# Install FEconv
USER root
WORKDIR /usr/local
RUN apt-get update --yes && \ 
  apt-get install -yq --no-install-recommends make wget zip unzip 
RUN wget https://github.com/victorsndvg/FEconv/archive/refs/heads/master.zip
RUN unzip master.zip && \
  rm master.zip
WORKDIR /usr/local/FEconv-master
RUN make -f Makefile.gfortran.linux
RUN ln -s /usr/local/FEconv-master/feconv /usr/local/bin/feconv

# Install locales
RUN apt-get update --yes && \ 
  apt-get install -yq --no-install-recommends locales-all && \
  locale

# cleanup
RUN apt-get --yes clean          && \
    apt-get --yes autoremove     && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* 

# Configure environment
ENV SHELL=/bin/bash \
    NB_USER="${NB_USER}" \
    NB_UID=${NB_UID} \
    NB_GID=${NB_GID} \
    LC_ALL=en_US.UTF-8 \
    LANG=en_US.UTF-8 \
    LANGUAGE=en_US.UTF-8 \
    HOME="/home/${NB_USER}" \
    HOSTNAME=jupyter

# Make sure the contents of our repo are in ${HOME}
COPY . ${HOME}
USER root
RUN chown -R ${NB_UID} ${HOME}

# Switch back to jovyan to avoid accidental container runs as root
USER ${NB_UID}

WORKDIR "${HOME}"