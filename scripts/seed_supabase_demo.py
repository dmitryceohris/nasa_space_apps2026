"""Publish the verified challenge snapshots to the same cache used by the web API."""
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from api import jobs, params, postgres, resolve


def seed():
    with postgres.connect() as conn:
        for path in sorted((ROOT / "frontend/challenge-demo").glob("*/api_response.json")):
            result = json.loads(path.read_text())
            images = []
            frames = result["frames"]
            query = result["query"]
            is_field = query["target"] is None
            p = params.parse(target=None if is_field else path.parent.name,
                             ra=query["ra"] if is_field else None, dec=query["dec"] if is_field else None,
                             det=query["det"], fov=query["fov"], mode=query["mode"], max_per_epoch=1 if is_field else 6)
            key = jobs._view_key(p, resolve.for_params(p, epoch=result.get("center_mjd")), source="archive")
            for index, frame in enumerate(frames):
                for field, suffix in (("image_url", ""), ("flags_url", "_flags")):
                    source = path.parent / f"frame_{index:02}{suffix}.png"
                    name = Path(frame[field]).name
                    frame[field] = f"/cache/{key}/{name}"
                    images.append((name, source.read_bytes()))
            postgres.save(conn, key, result, images, ttl=0)
            # Replays are explicitly labelled in their provenance and tied to a fixed epoch.
            conn.execute("update spherex.results set expires_at='infinity' where key=%s", (key,))
            print(path.parent.name, len(frames))


if __name__ == "__main__":
    seed()
