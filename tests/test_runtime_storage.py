"""Scratch quotas, job isolation, durable publication and scientific parity, without network."""
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
import errno
import hashlib
import io
import json
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits
from fastapi.testclient import TestClient

from api import cache, irsa, jobs, postgres
from api.main import app
from api.params import Params
from Exposure import calibration
from runtime_storage import StorageFull, current_workspace, worker_slot
from test_processing import CENTER, fake_exposure


@pytest.fixture
def scratch(monkeypatch, tmp_path):
    root = tmp_path / "cache"
    monkeypatch.setenv("SPHEREX_EPHEMERAL_STORAGE", "1")
    monkeypatch.setenv("CACHE_DIR", str(root))
    monkeypatch.setenv("SPHEREX_TMP_ROOT", str(tmp_path))
    monkeypatch.setenv("SPHEREX_CALIBRATION_CACHE", str(tmp_path / "calibration"))
    monkeypatch.setenv("SPHEREX_TMP_MIN_FREE_BYTES", "0")
    for name in ("XDG_CACHE_HOME", "XDG_CONFIG_HOME", "MPLCONFIGDIR"):
        monkeypatch.delenv(name, raising=False)
    return root


def test_disabled_mode_leaves_local_archive_untouched(monkeypatch, tmp_path):
    monkeypatch.delenv("SPHEREX_EPHEMERAL_STORAGE", raising=False)
    monkeypatch.setenv("CACHE_DIR", str(tmp_path / "missing"))
    with worker_slot() as workspace:
        assert workspace is None
    assert not (tmp_path / "missing").exists()


def test_platform_xdg_root_and_runtime_dependencies_are_never_evicted(scratch, monkeypatch):
    # Model Vercel's XDG_CACHE_HOME=/tmp with dependencies next to our private namespace.
    namespace = scratch.parent / "private"
    root = namespace / "cache"
    root.mkdir(parents=True)
    runtime = scratch.parent / "_vc_deps"
    runtime.mkdir()
    certificate = runtime / "certifi.pem"
    certificate.write_bytes(b"runtime certificate" * 100)
    monkeypatch.setenv("CACHE_DIR", str(root))
    monkeypatch.setenv("SPHEREX_TMP_ROOT", str(namespace))
    monkeypatch.setenv("XDG_CACHE_HOME", str(scratch.parent))
    monkeypatch.setenv("SPHEREX_TMP_MAX_BYTES", "100")
    with worker_slot() as workspace:
        workspace.atomic_write(workspace.path / "raw", [b"image"])
        assert workspace.roots == [root]
    assert certificate.read_bytes() == b"runtime certificate" * 100


def test_worker_slot_queues_competitors_and_cleans_dead_workspaces(scratch):
    orphan = scratch / "work/job-interrupted"
    orphan.mkdir(parents=True)
    (orphan / "partial.fits").write_bytes(b"interrupted")
    unrelated = scratch.parent / "personal.fits"
    unrelated.write_bytes(b"keep")
    with worker_slot() as workspace:
        assert not orphan.exists()
        assert current_workspace() is workspace
        assert workspace.path.is_dir()
        with worker_slot() as second:
            assert second is None
        active = workspace.path
    assert not active.exists() and current_workspace() is None
    assert unrelated.read_bytes() == b"keep"


def test_exception_and_timeout_remove_owned_files(scratch):
    with pytest.raises(TimeoutError):
        with worker_slot() as workspace:
            path = workspace.path
            (path / "raw.fits").write_bytes(b"partial")
            raise TimeoutError()
    assert not path.exists()
    with worker_slot() as next_workspace:
        assert next_workspace is not None


def test_quota_evicts_unowned_files_but_preserves_active_readers(scratch, monkeypatch):
    monkeypatch.setenv("SPHEREX_TMP_MAX_BYTES", "100")
    scratch.mkdir()
    old = scratch / "old.fits"
    old.write_bytes(b"x" * 80)
    with worker_slot() as workspace:
        with workspace.pin(old):
            with workspace.pin(old):
                pass
            with pytest.raises(StorageFull):
                workspace.atomic_write(workspace.path / "new.fits", [b"y" * 30])
            assert old.exists()
        workspace.atomic_write(workspace.path / "new.fits", [b"y" * 30])
        assert not old.exists()
        assert workspace.peak <= 100


def test_disk_full_and_checksum_failures_do_not_leave_partial_files(scratch):
    def interrupted():
        yield b"first chunk"
        raise OSError(errno.ENOSPC, "controlled full disk")
    with worker_slot() as workspace:
        dest = workspace.path / "raw.fits"
        with pytest.raises(StorageFull):
            workspace.atomic_write(dest, interrupted())
        with pytest.raises(ValueError, match="checksum"):
            workspace.atomic_write(dest, [b"bad"], expected_sha256="0" * 64)
        assert list(workspace.path.iterdir()) == []


def test_free_space_reserve_is_checked(scratch, monkeypatch):
    from runtime_storage import shutil
    monkeypatch.setenv("SPHEREX_TMP_MIN_FREE_BYTES", "100")
    monkeypatch.setattr(shutil, "disk_usage", lambda root: type("Usage", (), {"free":99})())
    with pytest.raises(StorageFull):
        with worker_slot():
            pytest.fail("started without free space reserve")


class Connection:
    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False

    def execute(self, *args):
        return self

    def fetchone(self):
        return (True,)

    def commit(self):
        pass


@pytest.mark.parametrize("failure", [None, "save", "commit"])
def test_postgres_keeps_images_until_commit_and_always_cleans(scratch, monkeypatch, failure):
    conn = Connection()
    state = {}
    monkeypatch.setattr(postgres, "connect", lambda: conn)
    monkeypatch.setattr(postgres, "read", lambda *args: None)
    def runner(job):
        workspace = current_workspace()
        state["path"] = workspace.path
        render = workspace.path / "render"
        render.mkdir()
        (render / "f0.png").write_bytes(b"png")
        return {"status":"ready", "frames":[{"image_url":"/cache/key/f0.png"}]}
    def save(conn, key, result, images, ttl):
        assert images == [("f0.png", b"png")]
        assert state["path"].exists()
        if failure == "save":
            raise RuntimeError("save failed")
    def commit():
        assert state["path"].exists()
        if failure == "commit":
            raise RuntimeError("commit failed")
    conn.commit = commit
    monkeypatch.setattr(postgres, "save", save)
    if failure:
        with pytest.raises(RuntimeError, match="failed"):
            postgres.run(jobs.Job("key", "id"), runner)
    else:
        assert postgres.run(jobs.Job("key", "id"), runner)[0] == 200
    assert not state["path"].exists()


def test_postgres_queue_does_not_block_cached_results(scratch, monkeypatch):
    monkeypatch.setattr(postgres, "connect", Connection)
    monkeypatch.setattr(postgres, "read", lambda *args: None)
    with worker_slot():
        status, payload = postgres.run(jobs.Job("key", "id"), lambda job: pytest.fail("busy runner started"))
        assert status == 202 and payload["stage"] == "queued"
        monkeypatch.setattr(postgres, "read", lambda *args: {"status":"ready"})
        assert postgres.run(jobs.Job("key", "id"), lambda job: pytest.fail("cached runner started"))[0] == 200


def test_storage_failure_is_public_503_and_not_a_ready_result(monkeypatch):
    def fail(*args):
        raise StorageFull()
    monkeypatch.setattr(jobs, "view", fail)
    with TestClient(app) as client:
        result = client.get("/api/view?target=barnard")
    assert result.status_code == 503
    assert result.json()["error"]["code"] == "storage_full"
    assert result.headers["retry-after"] == "5"


def test_pipeline_does_not_disguise_disk_full_as_missing_pixels(monkeypatch):
    monkeypatch.setattr(irsa, "exposures", lambda *args:[{"uri":"u0","mjd":60900.}])
    def full(*args, **kwargs):
        raise OSError(errno.ENOSPC, "controlled full disk")
    monkeypatch.setattr(irsa, "cutout", full)
    p = Params(None, CENTER.ra.deg, CENTER.dec.deg, 3., 3, "coadd", 6)
    r = {"center":(CENTER.ra.deg,CENTER.dec.deg),"motion":None,"label":"test"}
    with pytest.raises(StorageFull):
        jobs._pipeline(jobs.Job("key", "id"), p, r)
    monkeypatch.setattr(irsa, "cutout", lambda *args:"raw")
    monkeypatch.setattr(jobs, "load_exposure", full)
    with pytest.raises(StorageFull):
        jobs._pipeline(jobs.Job("key", "id"), p, r)


def test_spectrum_batches_preserve_the_same_points_and_counts(scratch, monkeypatch):
    from astropy.table import Table
    from processing import pipeline
    exps = [fake_exposure(CENTER,61000.+40*i,seed=i) for i in range(8)]
    rows = Table({"detector":[3]*8,"mjd":[e.mjd for e in exps]})
    monkeypatch.setattr(pipeline, "query_spherex", lambda *args:rows)
    batch_sizes = []
    def download(selected, *args, workspace=None, **kwargs):
        batch_sizes.append(len(selected))
        return [e for e in exps if e.mjd in [float(row["mjd"]) for row in selected]]
    monkeypatch.setattr(pipeline, "download_rows", download)
    baseline = pipeline.build_spectrum(CENTER, lambda mjd:CENTER)
    with worker_slot() as workspace:
        result = pipeline.build_spectrum(CENTER, lambda mjd:CENTER, workspace=workspace)
    assert result == baseline
    assert batch_sizes == [8,6,2]


@pytest.mark.parametrize("mode", ["single", "coadd", "diff", "trail"])
def test_ephemeral_frames_are_identical_and_release_each_epoch(scratch, monkeypatch, mode):
    monkeypatch.setattr(cache, "CACHE", scratch)
    exposures = {f"u{i}":fake_exposure(CENTER, mjd, seed=i)
                 for i,mjd in enumerate([60900.,60901.,61100.,61101.])}
    meta = [{"uri":name,"mjd":e.mjd} for name,e in exposures.items()]
    monkeypatch.setattr(irsa, "exposures", lambda *args:meta)
    def cutout(uri, *args, workspace=None):
        if workspace is None:
            return uri
        # The next epoch must not inherit raw files from the previous epoch.
        if uri == "u2":
            assert not (workspace.path / "u0").exists()
            assert not (workspace.path / "u1").exists()
        path = workspace.path / uri
        workspace.atomic_write(path, [b"fits"])
        return path
    monkeypatch.setattr(irsa, "cutout", cutout)
    monkeypatch.setattr(jobs, "load_exposure", lambda path, **kw:exposures[Path(path).name])
    p = Params(None, CENTER.ra.deg, CENTER.dec.deg, 3., 3, mode, 6)
    r = {"center":(CENTER.ra.deg,CENTER.dec.deg),"motion":None,"label":"test"}
    scratch.mkdir()
    baseline = jobs._pipeline(jobs.Job("key", "id"), p, r)
    expected = {path.name:path.read_bytes() for path in (scratch / "key").glob("*.png")}
    with worker_slot() as workspace:
        result = jobs._pipeline(jobs.Job("key", "id"), p, r)
        actual = {path.name:path.read_bytes() for path in (workspace.path / "render").glob("*.png")}
        assert result == baseline
        assert actual == expected
        assert not list(workspace.path.glob("u*"))


def test_calibration_misses_are_single_flight_and_corruption_is_repaired(scratch, monkeypatch, tmp_path):
    def payload(hdus):
        stream = io.BytesIO()
        fits.HDUList(hdus).writeto(stream)
        return stream.getvalue()
    header = fits.Header({"DETECTOR":3,"BUNIT":"arcsec2"})
    contents = {
        "wave":payload([fits.PrimaryHDU(),fits.ImageHDU(np.ones((3,3)),name="CWAVE"),
                        fits.ImageHDU(np.ones((3,3)),name="CBAND")]),
        "solid_angle":payload([fits.PrimaryHDU(),fits.ImageHDU(np.ones((3,3)),header=header)]),
    }
    row = {"detector":3,"releases":["dr1"],"pipeline_versions":["7.0.5"]}
    for kind,content in contents.items():
        row.update({kind+"_file":kind+".fits",kind+"_url":"https://test/"+kind,
                    kind+"_sha256":hashlib.sha256(content).hexdigest()})
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps({"calibrations":[row]}))
    calls = []
    @contextmanager
    def stream(method, url, **kwargs):
        calls.append(url)
        class Response:
            def raise_for_status(self):
                pass
            def iter_bytes(self, *args):
                return iter([contents[url.rsplit("/",1)[1]]])
        yield Response()
    monkeypatch.setattr("runtime_storage.httpx.stream", stream)
    calibration._load_registered_calibration.cache_clear()
    with worker_slot() as workspace:
        def read():
            with workspace.activate():
                return calibration._registered_calibration(str(manifest),1,3,"dr1","7.0.5")
        with ThreadPoolExecutor(6) as pool:
            values = list(pool.map(lambda _:read(), range(6)))
        assert len(calls) == 2 and all(v is values[0] for v in values)
        assert not list((tmp_path / "calibration").glob("*.fits"))
    calibration._load_registered_calibration.cache_clear()
    # The ordinary local disk cache still repairs corrupted downloads.
    (tmp_path / "calibration/wave.fits").write_bytes(b"broken")
    (tmp_path / "calibration/solid_angle.fits").write_bytes(contents["solid_angle"])
    local = calibration._registered_calibration(str(manifest),1,3,"dr1","7.0.5")
    assert local.detector == 3 and len(calls) == 3
    assert np.array_equal(local.cwave,values[0].cwave)
    assert np.array_equal(local.cband,values[0].cband)
    assert np.array_equal(local.solid_angle,values[0].solid_angle)
    assert local.provenance == values[0].provenance
    assert not list((tmp_path / "calibration").glob(".part-*"))
    calibration._load_registered_calibration.cache_clear()
