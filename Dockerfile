FROM debian:8.5

MAINTAINER Karsten Wiesner <karsten.wiesner@aei.mpg.de>

# Setup environment vars
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV CONDA_DIR=/opt/conda \
    SHELL=/bin/bash \
    NB_USER=lisa \
    NB_UID=1000 \
    NB_GID=100

ENV PATH=$CONDA_DIR/bin:$PATH \
    HOME=/home/$NB_USER

# Add bash script from the host directory
ADD fix-permissions /usr/local/bin/fix-permissions

# Create user with UID=1000 and in the 'users' group
# and make sure these dirs are writable by the `users` group.
RUN useradd -m -s /bin/bash -N -u $NB_UID $NB_USER && \
    mkdir -p $CONDA_DIR && \
    chown $NB_USER:$NB_GID $CONDA_DIR && \
    fix-permissions $HOME && \
    fix-permissions $CONDA_DIR

# Configure apt-get
RUN apt-get update --fix-missing && apt-get install -y --no-install-recommends apt-utils

# Install packages
RUN apt-get install -y wget bzip2 ca-certificates \
    libglib2.0-0 libxext6 libsm6 libxrender1 \
    git python2.7 python2.7-dev python-pip swig sudo

# Configure python
RUN update-alternatives --install /usr/bin/python python /usr/bin/python2.7 20
RUN update-alternatives --set python /usr/bin/python2.7

# Install python modules
RUN pip install numpy
RUN pip install pyRXP

# Give NB_USER sudo rights (and reset user password to user-name)
RUN echo "$NB_USER:$NB_USER" | chpasswd
RUN echo "$NB_USER ALL=(ALL) ALL" |cat >> /etc/sudoers

# Setup tini
RUN apt-get install -y curl grep sed dpkg && \
    TINI_VERSION=`curl https://github.com/krallin/tini/releases/latest | grep -o "/v.*\"" | sed 's:^..\(.*\).$:\1:'` && \
    curl -L "https://github.com/krallin/tini/releases/download/v${TINI_VERSION}/tini_${TINI_VERSION}.deb" > tini.deb && \
    dpkg -i tini.deb && \
    rm tini.deb && \
    apt-get clean

ENTRYPOINT [ "/usr/bin/tini", "--" ]

USER $NB_USER
WORKDIR  $HOME
CMD [ "/bin/bash" ]
