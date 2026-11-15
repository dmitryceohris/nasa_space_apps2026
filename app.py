"""Vercel entry point for the SPHEREx viewer."""
import os
from pathlib import Path

os.environ.setdefault("CACHE_DIR", "/tmp/spherex")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("PREWARM", "0")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/spherex-xdg-cache")
os.environ.setdefault("XDG_CONFIG_HOME", "/tmp/spherex-xdg-config")
os.environ.setdefault("MAX_JOB_SECONDS", "260")

os.environ.setdefault("SPHEREX_CALIBRATION_MANIFEST",
                      str(Path(__file__).parent / "data/validation/science-remediation/calibrations.json"))
os.environ.setdefault("SPHEREX_CALIBRATION_CACHE", "/tmp/spherex-calibrations")

from api.main import app
