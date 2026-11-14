"""SPHEREx blink API."""
import sys
from pathlib import Path

# Exposure and processing live in src/. pytest adds it through pyproject.toml; uvicorn does not, so do it here.
_SRC = str(Path(__file__).resolve().parent.parent / "src")
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
