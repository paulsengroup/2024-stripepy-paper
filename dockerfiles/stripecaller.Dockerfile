# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM python:3.12 AS base

ARG CONTAINER_VERSION=0.1.0

RUN apt-get update \
&& apt-get install -q -y --no-install-recommends time \
&& rm -rf /var/lib/apt/lists/*

RUN python3 -m venv /opt/stripecaller --upgrade-deps \
&& /opt/stripecaller/bin/pip install --no-cache-dir \
    cooler \
    matplotlib \
    'numpy<2' \
    scipy \
    "stripecaller==$CONTAINER_VERSION"

# Populate bytecode cache
ENV PYTHONDONTWRITEBYTECODE=
ENV MPLCONFIGDIR=/tmp/.matplotlib
RUN mkdir "$MPLCONFIGDIR"

RUN /opt/stripecaller/bin/call-stripes --help

WORKDIR /data
ENTRYPOINT ["/opt/stripecaller/bin/call-stripes"]
ENV PATH="$PATH:/opt/stripecaller/bin"
