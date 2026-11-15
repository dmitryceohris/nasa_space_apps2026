"""Bounded scratch space for serverless jobs; ordinary local archives are untouched."""
from contextlib import contextmanager
from contextvars import ContextVar
from collections import Counter
import errno
import fcntl
import hashlib
import logging
import os
from pathlib import Path
import shutil
import tempfile
import threading
import time

import httpx

log = logging.getLogger(__name__)
_current = ContextVar("spherex_workspace", default=None)
STORAGE_MESSAGE = "Недостаточно временного места для обработки снимков. Повторите запрос позже."


class StorageFull(OSError):
    def __init__(self):
        super().__init__(errno.ENOSPC, STORAGE_MESSAGE)


def raise_if_storage_error(error):
    if isinstance(error, StorageFull):
        raise error
    if isinstance(error, OSError) and error.errno in (errno.ENOSPC, errno.EDQUOT):
        raise StorageFull() from error


def enabled():
    return os.environ.get("SPHEREX_EPHEMERAL_STORAGE") == "1"


def current_workspace():
    return _current.get()


class Workspace:
    def __init__(self, cache_root):
        self.root = cache_root
        namespace = Path(os.environ.get("SPHEREX_TMP_ROOT", cache_root)).resolve()
        if not cache_root.resolve().is_relative_to(namespace):
            raise ValueError("CACHE_DIR must belong to the application's scratch namespace")
        self.roots = list(dict.fromkeys([cache_root, *(
            Path(os.environ[name]) for name in
            ("SPHEREX_CALIBRATION_CACHE", "XDG_CACHE_HOME", "XDG_CONFIG_HOME", "MPLCONFIGDIR")
            if os.environ.get(name) and Path(os.environ[name]).resolve().is_relative_to(namespace)
            and Path(os.environ[name]).resolve() != namespace)]))
        self.maximum = int(os.environ.get("SPHEREX_TMP_MAX_BYTES", "350000000"))
        self.reserve = int(os.environ.get("SPHEREX_TMP_MIN_FREE_BYTES", "75000000"))
        self.guard = threading.RLock()
        self.pinned = Counter()
        self.peak = 0
        self.path = None

    def files(self):
        # Resolve duplicates where a configured cache lives inside another managed root.
        return {p for root in self.roots if root.is_dir() for p in root.rglob("*")
                if p.is_file() and not p.is_symlink()}

    def ensure_space(self, incoming=0):
        with self.guard:
            files = self.files()
            sizes, modified = {}, {}
            for path in files:
                try:
                    stat = path.stat()
                    sizes[path], modified[path] = stat.st_size, stat.st_mtime
                except FileNotFoundError:
                    continue
            used = sum(sizes.values())
            free = shutil.disk_usage(self.root).free
            candidates = sorted((p for p in sizes if p not in self.pinned
                                 and (self.path is None or not p.is_relative_to(self.path))),
                                key=modified.__getitem__)
            while used + incoming > self.maximum or free - incoming < self.reserve:
                if not candidates:
                    raise StorageFull()
                path = candidates.pop(0)
                path.unlink(missing_ok=True)
                used -= sizes[path]
                free = shutil.disk_usage(self.root).free
            self.peak = max(self.peak, used + incoming)

    @contextmanager
    def activate(self):
        token = _current.set(self)
        try:
            yield self
        finally:
            _current.reset(token)

    @contextmanager
    def pin(self, *paths):
        paths = tuple(map(Path, paths))
        with self.guard:
            self.pinned.update(paths)
        try:
            yield
        finally:
            with self.guard:
                self.pinned.subtract(paths)
                self.pinned += Counter()  # Drop zero counts after nested/concurrent readers finish.

    def atomic_write(self, destination, chunks, expected_sha256=None):
        destination = Path(destination)
        destination.parent.mkdir(parents=True, exist_ok=True)
        temporary = None
        try:
            with self.pin(destination):
                self.ensure_space()
                with tempfile.NamedTemporaryFile(dir=destination.parent, prefix=".part-", delete=False,
                                                 buffering=0) as stream:
                    temporary = Path(stream.name)
                    digest = hashlib.sha256()
                    with self.pin(temporary):
                        for chunk in chunks:
                            with self.guard:
                                self.ensure_space(len(chunk))
                                stream.write(chunk)
                            digest.update(chunk)
                        if expected_sha256 and digest.hexdigest() != expected_sha256:
                            raise ValueError("calibration_checksum_mismatch")
                        temporary.replace(destination)
        except OSError as error:
            raise_if_storage_error(error)
            raise
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)
        return destination

    def download(self, url, destination=None, *, expected_sha256=None, timeout=45):
        if destination is None:
            destination = self.path / (hashlib.sha256(url.encode()).hexdigest() + ".fits")
        for attempt in range(3):
            try:
                with httpx.stream("GET", url, timeout=timeout, follow_redirects=True) as response:
                    response.raise_for_status()
                    return self.atomic_write(destination, response.iter_bytes(128 * 1024), expected_sha256)
            except httpx.HTTPError as error:
                code = error.response.status_code if isinstance(error, httpx.HTTPStatusError) else None
                if attempt == 2 or (code is not None and code < 500 and code not in (408, 429)):
                    raise
                time.sleep(2 ** attempt)

    def verified_bytes(self, url, expected_sha256):
        """Calibration arrays are cached in RAM, avoiding another 50 MB per detector on disk."""
        for attempt in range(3):
            try:
                data, digest = bytearray(), hashlib.sha256()
                with httpx.stream("GET", url, timeout=45, follow_redirects=True) as response:
                    response.raise_for_status()
                    for chunk in response.iter_bytes(128 * 1024):
                        if len(data) + len(chunk) > 64_000_000:
                            raise ValueError("calibration_size_limit")
                        data.extend(chunk)
                        digest.update(chunk)
                if digest.hexdigest() != expected_sha256:
                    raise ValueError("calibration_checksum_mismatch")
                return bytes(data)
            except httpx.HTTPError as error:
                code = error.response.status_code if isinstance(error, httpx.HTTPStatusError) else None
                if attempt == 2 or (code is not None and code < 500 and code not in (408, 429)):
                    raise
                time.sleep(2 ** attempt)


@contextmanager
def worker_slot(*, blocking=False):
    """One heavy job per scratch filesystem, across threads and processes."""
    if not enabled():
        yield None
        return
    root = Path(os.environ.get("CACHE_DIR", "/tmp/spherex"))
    root.mkdir(parents=True, exist_ok=True)
    # Outside eviction roots: removing a locked file would allow a second lock inode.
    lock_path = root.parent / ("." + root.name + "-worker.lock")
    with lock_path.open("a+b") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | (0 if blocking else fcntl.LOCK_NB))
        except BlockingIOError:
            yield None
            return
        workspace = Workspace(root)
        try:
            log.info("scratch job starting: namespace=%s free_bytes=%d", root, shutil.disk_usage(root).free)
            jobs = root / "work"
            jobs.mkdir(exist_ok=True)
            # All workspaces share this lock; these can only belong to dead invocations.
            for orphan in jobs.glob("job-*"):
                if orphan.is_dir() and not orphan.is_symlink():
                    shutil.rmtree(orphan)
            workspace.ensure_space()
            workspace.path = Path(tempfile.mkdtemp(prefix="job-", dir=jobs))
            with workspace.activate():
                yield workspace
        except OSError as error:
            raise_if_storage_error(error)
            raise
        finally:
            if workspace.path is not None:
                shutil.rmtree(workspace.path)
            log.info("scratch job finished: peak_managed_bytes=%d limit_bytes=%d free_bytes=%d",
                     workspace.peak, workspace.maximum, shutil.disk_usage(root).free)
            fcntl.flock(lock, fcntl.LOCK_UN)
