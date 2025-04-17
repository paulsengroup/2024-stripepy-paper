# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

ARG CONTAINER_VERSION=1.0.0

FROM paulsengroup/stripepy:${CONTAINER_VERSION} AS base

RUN apt-get update \
&& apt-get install -q -y --no-install-recommends time \
&& rm -rf /var/lib/apt/lists/*

RUN stripepy call --help
RUN stripepy --version

WORKDIR /data
ENV PYTHONDONTWRITEBYTECODE=1
