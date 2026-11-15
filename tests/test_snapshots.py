"""The five saved targets must support every mode without a database or archive."""
import io

import pytest
from fastapi.testclient import TestClient
from PIL import Image

from api import jobs, resolve
from api.main import app


@pytest.fixture
def client(monkeypatch):
    monkeypatch.delenv("SPHEREX_DATABASE_URL", raising=False)
    monkeypatch.setattr(resolve, "_fetch_sesame", lambda *a: pytest.fail("snapshot used network"))
    monkeypatch.setattr(jobs, "_pipeline", lambda *a: pytest.fail("snapshot reprocessed FITS at request time"))
    with TestClient(app) as client:
        yield client


@pytest.mark.parametrize("target", resolve.PRESETS, ids=lambda t: t["id"])
@pytest.mark.parametrize("mode", ["single", "coadd", "diff", "trail"])
@pytest.mark.parametrize("maximum", [6, 1])
def test_saved_targets_have_real_mode_specific_frames(client, target, mode, maximum):
    query = {"target": target["id"], "det": target["det"], "fov": target["fov"],
             "mode": mode, "snapshot": "2026-11-15", "max_per_epoch": maximum}
    response = client.get("/api/view", params=query)
    assert response.status_code == 200, response.text
    result = response.json()
    assert result["query"]["mode"] == mode
    assert result["archive_source"] == "saved_sia2_replay"
    assert len(result["frames"]) >= 2
    assert result["exposure_counts"]["failed"] == 0
    if mode == "single":
        assert all(frame["n_exp"] == 1 for frame in result["frames"])
    if mode == "diff":
        assert result["processing"]["difference"]["name"] == "leave_one_out_median"
        assert all(frame["diff_noise_scale"] is not None for frame in result["frames"])
    for frame in result["frames"]:
        for field in ["image_url", "flags_url"]:
            image = client.get(frame[field])
            assert image.status_code == 200
            assert image.content.startswith(b"\x89PNG")
    if mode in {"diff", "trail"}:
        trail = client.get(result["trail_url"])
        assert trail.status_code == 200 and trail.content.startswith(b"\x89PNG")
    gif = client.get("/api/export.gif", params=query)
    assert gif.status_code == 200
    with Image.open(io.BytesIO(gif.content)) as image:
        assert image.n_frames == len(result["frames"])


def test_snapshot_does_not_substitute_another_field(client):
    response = client.get("/api/view", params={"target": "barnard", "det": 6, "fov": 20,
                                              "snapshot": "2026-11-15"})
    assert response.status_code == 404
    assert response.json()["error"]["code"] == "no_snapshot"


def test_unknown_snapshot_date_is_rejected(client):
    response = client.get("/api/view", params={"target": "barnard", "snapshot": "2026-01-01"})
    assert response.status_code == 400


@pytest.mark.parametrize("mode", ["single", "coadd", "diff", "trail"])
def test_exact_coordinate_snapshot_keeps_fixed_field_semantics(client, mode):
    response = client.get("/api/view", params={"ra": "08h55m10.831s", "dec": "-7d14m42.53s",
        "det": 6, "fov": 6, "mode": mode, "max_per_epoch": 1, "snapshot": "2026-11-15"})
    assert response.status_code == 200
    result = response.json()
    assert result["target_info"] is None
    assert all(f["predicted_px"] is None for f in result["frames"])


def test_packaged_snapshot_and_gif_do_not_open_configured_database(client, monkeypatch):
    from api import postgres
    monkeypatch.setenv("SPHEREX_DATABASE_URL", "unused")
    monkeypatch.setattr(postgres, "connect", lambda: pytest.fail("packaged snapshot used database"))
    query = {"target": "barnard", "det": 3, "fov": 3, "mode": "trail", "snapshot": "2026-11-15"}
    assert client.get("/api/view", params=query).status_code == 200
    assert client.get("/api/export.gif", params=query).status_code == 200


def test_incompatible_snapshot_processing_version_is_not_served(client, monkeypatch):
    from api import cache
    monkeypatch.setattr(cache, "VERSION", 999)
    response = client.get("/api/view", params={"target": "barnard", "det": 3, "fov": 3,
                                              "mode": "single", "snapshot": "2026-11-15"})
    assert response.status_code == 404
