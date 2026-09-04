FROM condaforge/miniforge3:24.11.3-0

COPY envs/pacusage.yml /tmp/pacusage.yml
RUN /opt/conda/bin/conda env create -n pacusage -f /tmp/pacusage.yml \
    && /opt/conda/bin/conda clean --all --yes

COPY pyproject.toml README.md LICENSE /opt/pacusage/
COPY src /opt/pacusage/src
RUN /opt/conda/envs/pacusage/bin/pip install \
    --no-deps --no-build-isolation /opt/pacusage

ENV PATH=/opt/conda/envs/pacusage/bin:$PATH
ENV LC_ALL=C
ENTRYPOINT ["/usr/bin/tini", "--"]
