FROM ghcr.io/prefix-dev/pixi:latest

WORKDIR /app
COPY pixi.toml .
RUN pixi install && pixi clean cache -y

ENV PATH=/app/.pixi/envs/default/bin:$PATH
