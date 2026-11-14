#!/bin/bash
# Demo mode: fixed query-centre date, no startup prewarm, ready-made cache.
# Usage: bash start_demo.sh [CACHE_DIR] [CENTER_MJD]
# The cache is only valid for the CENTER_MJD it was built with (the cache key depends on the propagated centre).
D="$(cd "$(dirname "$0")" && pwd)"
export CACHE_DIR="${1:-$D/cache}"
export CENTER_MJD="${2:-${CENTER_MJD:?set CENTER_MJD or pass it as the second argument}}"
export PREWARM=0
exec "$D/.venv/bin/uvicorn" --app-dir "$D" api.main:app --host "${HOST:-127.0.0.1}" --port "${PORT:-8000}"
