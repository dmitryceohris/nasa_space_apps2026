# Build: docker build -t spherex-api .
# Run:   docker run --rm -p 8000:8000 -v spherex-cache:/data/cache spherex-api
FROM python:3.12.13-slim

COPY --from=ghcr.io/astral-sh/uv:0.9.5 /uv /usr/local/bin/uv

WORKDIR /app
ENV UV_COMPILE_BYTECODE=1 UV_LINK_MODE=copy PYTHONUNBUFFERED=1 \
    CACHE_DIR=/data/cache HOST=0.0.0.0 PORT=8000

COPY pyproject.toml uv.lock .python-version ./
RUN uv sync --locked --no-dev --no-install-project

COPY . .

VOLUME /data/cache
EXPOSE 8000
HEALTHCHECK --interval=30s --timeout=5s --start-period=20s \
    CMD python -c "import urllib.request,sys; sys.exit(urllib.request.urlopen('http://127.0.0.1:8000/api/health').status != 200)"
CMD ["sh", "-c", "exec /app/.venv/bin/uvicorn --app-dir /app api.main:app --host $HOST --port $PORT"]
