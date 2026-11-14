#!/bin/bash
# Demo mode: fixed query-centre date, no startup prewarm, ready-made cache.
# Usage: bash start_demo.sh [CACHE_DIR] [CENTER_MJD]
# The cache is only valid for the CENTER_MJD it was built with (the cache key depends on the propagated centre).
set -e
D="$(cd "$(dirname "$0")" && pwd)"
if [ ! -x "$D/.venv/bin/uvicorn" ]; then
    uv sync --project "$D" --locked
fi
export CACHE_DIR="${1:-${CACHE_DIR:-$D/cache}}"
export CENTER_MJD="${2:-${CENTER_MJD:-61358}}"
export PREWARM=0
exec "$D/.venv/bin/uvicorn" --app-dir "$D" api.main:app --host "${HOST:-127.0.0.1}" --port "${PORT:-8000}"
