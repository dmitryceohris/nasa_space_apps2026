#!/bin/bash
# Set PREWARM=0 to skip warming the 5 presets at startup; CENTER_MJD pins the date the query centre is propagated to.
D="$(cd "$(dirname "$0")" && pwd)"
exec "$D/.venv/bin/uvicorn" --app-dir "$D" api.main:app --host "${HOST:-127.0.0.1}" --port "${PORT:-8000}"
