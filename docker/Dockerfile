ARG OS_VERSION="20.04"
ARG CUDA_VERSION="11.1.1"
FROM nvidia/cuda:${CUDA_VERSION}-devel-ubuntu${OS_VERSION}

USER claire-user

WORKDIR /home/claire-user

# Install package dependencies
RUN apt-get update -qq \
      && apt-get -y --no-install-recommends install \
         build-essential \
         ca-certificates \
         gdb \
         gfortran \
         wget \
         vim \
         cmake \
         python2 \
         zlib1g \
         zlib1g-dev \
         git \
         openmpi-bin \
         openmpi-doc \
         libopenmpi-dev
      && apt-get clean all \
      && rm -r /var/lib/apt/lists/*

RUN cd /home/claire-user && git clone https://github.com/andreasmang/claire.git && cd claire && git checkout single_gpu

# Install dependencies
RUN cd /user/claire-user/claire/deps/ && make && source env_source.sh 

# Build claire
RUN cd /user/claire-user/claire/ && make

CMD ["/bin/bash"]
