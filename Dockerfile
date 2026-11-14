FROM python:3.12.13-slim
WORKDIR /app
RUN pip install --no-cache-dir uv==0.11.29
COPY pyproject.toml uv.lock ./
RUN uv sync --locked --no-dev --no-install-project
COPY . .
RUN useradd --uid 10001 --create-home app && mkdir /cache && chown app:app /cache
ENV PATH="/app/.venv/bin:$PATH" \
    CACHE_DIR=/cache PREWARM=0 HOST=0.0.0.0 PORT=8000 \
    PYTHONUNBUFFERED=1 MPLCONFIGDIR=/tmp/matplotlib
USER app
EXPOSE 8000
VOLUME ["/cache"]
HEALTHCHECK --interval=30s --timeout=5s --start-period=20s \
    CMD python -c "import os,urllib.request; urllib.request.urlopen('http://127.0.0.1:'+os.environ.get('PORT','8000')+'/api/health', timeout=3)"
CMD ["bash", "start.sh"]
