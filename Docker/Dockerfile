FROM ubuntu:18.04

RUN apt-get update && apt-get install -y \
    zsh \
    git \
    cmake \
    clang \
    libboost-all-dev \
    g++ \
    wget

# getting Gudhi 1.3
RUN wget https://gforge.inria.fr/frs/download.php/file/37129/2016-09-12-16-12-33_GUDHI_1.3.1.tar.gz
RUN tar -xzf 2016-09-12-16-12-33_GUDHI_1.3.1.tar.gz && \
    cd 2016-09-12-16-12-33_GUDHI_1.3.1 && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && \
    make install

# getting eigen 3.1.2
RUN wget https://github.com/eigenteam/eigen-git-mirror/archive/3.1.2.tar.gz && \
    tar -xzf /3.1.2.tar.gz && \
    cd eigen-git-mirror-3.1.2 && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && \
    make install

RUN git clone https://github.com/crvs/coeff-flow && \
    cd coeff-flow && \
    git submodule update --init --recursive && \
    mkdir build && \
    cd build && \
    cmake -DCMAKE_MODULE_PATH="/eigen-git-mirror-3.1.2/cmake" ..

# RUN cd coeff-flow/build && \
    # make

CMD ["zsh"]
