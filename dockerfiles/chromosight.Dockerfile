# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM python:3.12 AS base

ARG CONTAINER_VERSION=1.6.3

RUN apt-get update \
&& apt-get install -q -y --no-install-recommends time \
&& rm -rf /var/lib/apt/lists/*

RUN python3 -m venv /opt/chromosight --upgrade-deps \
&& /opt/chromosight/bin/pip install "chromosight==$CONTAINER_VERSION" 'numpy<2' --no-cache-dir

# Populate bytecode cache
ENV PYTHONDONTWRITEBYTECODE=

RUN /opt/chromosight/bin/chromosight detect --help
RUN /opt/chromosight/bin/chromosight --version

WORKDIR /data
ENTRYPOINT ["/opt/chromosight/bin/chromosight"]
ENV PATH="$PATH:/opt/chromosight/bin"
