import sys
from pathlib import Path

# src/ holds the shared libraries (Exposure: FITS loading, processing: trail/photometry). Make them importable
# however the server is started (uvicorn from the repo root, start.sh, tests).
_SRC = str(Path(__file__).resolve().parent.parent / "src")
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
