# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM python:3.11 AS base

ARG CONTAINER_VERSION=1.1.65.18

RUN apt-get update \
&& apt-get install -q -y --no-install-recommends libgl1-mesa-glx time \
&& rm -rf /var/lib/apt/lists/*

RUN python3 -m venv /opt/stripenn --upgrade-deps \
&& /opt/stripenn/bin/pip install "stripenn==$CONTAINER_VERSION" 'numpy<2' --no-cache-dir

# Populate bytecode cache
ENV PYTHONDONTWRITEBYTECODE=
ENV MPLCONFIGDIR=/tmp/.matplotlib
RUN mkdir "$MPLCONFIGDIR"

RUN /opt/stripenn/bin/stripenn compute --help
RUN /opt/stripenn/bin/python -c 'from importlib.metadata import version; print(version("stripenn"))'

WORKDIR /data
ENTRYPOINT ["/opt/stripenn/bin/stripenn"]
ENV PATH="$PATH:/opt/stripenn/bin"
