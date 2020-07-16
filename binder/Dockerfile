FROM atlasamglab/stats-base:root6.22.00-python3.8

# Run the following commands as super user (root):
USER root
WORKDIR /root

SHELL [ "/bin/bash", "-c" ]

# Install required packages for notebooks
RUN apt-get update && \
    apt-get install -qq -y \
        wget curl \
        git && \
    apt-get -y autoclean && \
    apt-get -y autoremove && \
    rm -rf /var/lib/apt/lists/* && \
    rm -rf /tmp/*

# c.f. https://mybinder.readthedocs.io/en/latest/tutorials/dockerfile.html#preparing-your-dockerfile
ARG NB_USER=jovyan
ARG NB_UID=1000
ENV USER ${NB_USER}
ENV NB_UID ${NB_UID}
ENV HOME /home/${NB_USER}

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}

RUN python -m pip install --upgrade --no-cache-dir pip setuptools wheel && \
    python -m pip install --no-cache-dir \
       jupyter \
       metakernel \
       zmq

# Make sure the contents of the repo are in ${HOME}
COPY . ${HOME}
USER root
RUN python -m pip install --no-cache-dir -r ${HOME}/requirements.txt
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}
WORKDIR ${HOME}

# Specify the default command to run
ENTRYPOINT []
CMD ["jupyter", "notebook", "--ip", "0.0.0.0"]
