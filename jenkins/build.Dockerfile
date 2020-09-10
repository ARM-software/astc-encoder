FROM ubuntu:18.04

LABEL build.environment.version="2.0.0"

RUN useradd -u 1001 -U -m -c Jenkins jenkins

RUN apt update && apt-get install -y \
    clang \
    g++ \
    gcc \
    git \
    make \
    python3 \
    python3-junit.xml \
    python3-numpy \
    python3-pil
