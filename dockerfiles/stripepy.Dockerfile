# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

ARG CONTAINER_VERSION=1.0.0

FROM paulsengroup/stripepy:${CONTAINER_VERSION} AS base

RUN apt-get update \
&& apt-get install -q -y --no-install-recommends time \
&& rm -rf /var/lib/apt/lists/*

# Populate bytecode cache
ENV PYTHONDONTWRITEBYTECODE=
ENV MPLCONFIGDIR=/tmp/.matplotlib
RUN mkdir "$MPLCONFIGDIR"

RUN stripepy --help \
&& stripepy call --help \
&& stripepy --version

WORKDIR /data
