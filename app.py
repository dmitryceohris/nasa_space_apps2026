"""Vercel entry point for the SPHEREx viewer."""
import os
from pathlib import Path

# Vercel can predefine XDG_CACHE_HOME=/tmp, alongside its hydrated Python runtime.
# Never treat that shared directory as an application cache or an eviction root.
_scratch = Path("/tmp/spherex-runtime")
os.environ["SPHEREX_TMP_ROOT"] = str(_scratch)
os.environ["CACHE_DIR"] = str(_scratch / "cache")
os.environ["MPLCONFIGDIR"] = str(_scratch / "matplotlib")
os.environ.setdefault("PREWARM", "0")
os.environ["XDG_CACHE_HOME"] = str(_scratch / "libraries/cache")
os.environ["XDG_CONFIG_HOME"] = str(_scratch / "libraries/config")
os.environ.setdefault("MAX_JOB_SECONDS", "260")
os.environ.setdefault("SPHEREX_EPHEMERAL_STORAGE", "1")

os.environ.setdefault("SPHEREX_CALIBRATION_MANIFEST",
                      str(Path(__file__).parent / "data/validation/science-remediation/calibrations.json"))
os.environ["SPHEREX_CALIBRATION_CACHE"] = str(_scratch / "calibrations")

from api.main import app
