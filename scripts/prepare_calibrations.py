"""Download the pinned DR1 calibration baseline; preserve hashes and source URLs.

These are the six products per type listed by IRSA's public S3 inventory.
No latest-by-date selection or implicit release mixing is performed.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path

import httpx

BASE = "https://nasa-irsa-spherex.s3.us-east-1.amazonaws.com/"


def prepare(directory):
    directory.mkdir(parents=True,exist_ok=True)
    def detector(n):
        row = {"detector":n,"releases":["dr1","qr3"],"pipeline_versions":["7.0.5"]}
        for name,kind,version in [("wave","spectral_wcs","cal-swcs-v5-2026-191"),
                                  ("solid_angle","solid_angle_pixel_map","cal-sapm-v3-2026-191")]:
            filename = f"{kind}_D{n}_spx_{version}.fits"
            url = f"{BASE}dr1/{kind}/{version}/{n}/{filename}"
            path = directory / filename
            if not path.exists():
                with httpx.stream("GET",url,timeout=60,follow_redirects=True) as response:
                    response.raise_for_status()
                    with path.with_suffix(".part").open("wb") as out:
                        for data in response.iter_bytes(): out.write(data)
                path.with_suffix(".part").replace(path)
            row.update({f"{name}_file":filename,f"{name}_url":url,
                        f"{name}_sha256":hashlib.sha256(path.read_bytes()).hexdigest()})
        print(f"D{n}: downloaded and hashed",flush=True)
        return row
    with ThreadPoolExecutor(max_workers=3) as pool:
        rows = list(pool.map(detector,range(1,7)))
    result = {"source":"IRSA DR1 release calibration inventory, accessed 2026-11-15",
              "matching_policy":"explicit release and pipeline version; never inferred from observation date",
              "calibrations":rows}
    (directory / "manifest.json").write_text(json.dumps(result,indent=2)+"\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=Path("cache/calibrations"))
    prepare(parser.parse_args().output)
