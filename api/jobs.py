"""Фоновые задачи с кэшем на диске: «получить готовое или запустить и отдать прогресс».

Скачивание с IRSA занимает минуты, поэтому HTTP-запрос не ждёт: первый вызов запускает
задачу и отвечает 202 с прогрессом, фронтенд повторяет тот же запрос, пока не получит 200.
Готовый результат лежит в CACHE_DIR/<kind>/<key>/result.json и переживает перезапуск сервера.
"""
from __future__ import annotations

import hashlib
import json
import logging
import os
import threading
import traceback
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, dataclass, field
from pathlib import Path

log = logging.getLogger("spherex.jobs")

CACHE_DIR = Path(os.environ.get("CACHE_DIR", Path(__file__).resolve().parents[1] / "cache"))


def cache_key(**params) -> str:
    raw = json.dumps(params, sort_keys=True, default=str)
    return hashlib.sha1(raw.encode()).hexdigest()[:16]


@dataclass
class JobState:
    status: str = "running"            # running | error
    stage: str = "queued"              # queued | query | download | process | render
    done: int = 0
    total: int = 0
    message: str = "В очереди…"
    error: str | None = field(default=None, repr=False)


class JobManager:
    def __init__(self, workers: int = 2):
        self._pool = ThreadPoolExecutor(max_workers=workers, thread_name_prefix="job")
        self._jobs: dict[str, JobState] = {}
        self._lock = threading.Lock()

    def job_dir(self, kind: str, key: str) -> Path:
        return CACHE_DIR / kind / key

    def get_or_start(self, kind: str, key: str,
                     fn: Callable[[Path, Callable], dict]) -> tuple[int, dict]:
        """(200, результат) если готово, иначе (202, прогресс). Ошибка -> (200, status=error) один раз."""
        d = self.job_dir(kind, key)
        result = d / "result.json"
        if result.exists():
            return 200, json.loads(result.read_text(encoding="utf-8"))
        jid = f"{kind}/{key}"
        with self._lock:
            st = self._jobs.get(jid)
            if st is not None and st.status == "error":
                self._jobs.pop(jid)  # отдать ошибку один раз; следующий запрос перезапустит
                return 200, {"status": "error", "message": st.error or st.message}
            if st is None:
                st = JobState()
                self._jobs[jid] = st
                self._pool.submit(self._run, jid, d, fn, st)
        return 202, {"status": "running", "job_id": key, **{k: v for k, v in asdict(st).items()
                                                             if k in ("stage", "done", "total", "message")}}

    def _run(self, jid: str, d: Path, fn, st: JobState) -> None:
        def progress(stage: str, done: int, total: int, message: str) -> None:
            st.stage, st.done, st.total, st.message = stage, done, total, message

        try:
            d.mkdir(parents=True, exist_ok=True)
            res = fn(d, progress)
            tmp = d / "result.json.tmp"
            tmp.write_text(json.dumps(res, ensure_ascii=False), encoding="utf-8")
            tmp.replace(d / "result.json")  # атомарно: фронт не увидит недописанный файл
            with self._lock:
                self._jobs.pop(jid, None)
            log.info("job %s ready", jid)
        except Exception as e:  # noqa: BLE001
            log.error("job %s failed: %s\n%s", jid, e, traceback.format_exc())
            st.status, st.error = "error", f"{type(e).__name__}: {e}"

    def wait_all(self, timeout: float = 600) -> None:
        """Для тестов: дождаться всех задач."""
        import time
        t0 = time.time()
        while time.time() - t0 < timeout:
            with self._lock:
                if not any(s.status == "running" for s in self._jobs.values()):
                    return
            time.sleep(0.1)
