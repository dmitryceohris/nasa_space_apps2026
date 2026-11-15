"""Shared results and rendered images for serverless deployments."""
import os
from pathlib import Path

import psycopg
from psycopg.types.json import Jsonb

from . import cache
from runtime_storage import enabled, worker_slot


def connect():
    return psycopg.connect(os.environ["SPHEREX_DATABASE_URL"], connect_timeout=8,
                           prepare_threshold=None)


def read(conn, key):
    row = conn.execute(
        "select payload from spherex.results where key=%s and version=%s and expires_at > now()",
        (key, cache.VERSION),
    ).fetchone()
    return row[0] if row else None


def read_image(key, filename):
    with connect() as conn:
        row = conn.execute("select png from spherex.images where result_key=%s and filename=%s",
                           (key, filename)).fetchone()
        return bytes(row[0]) if row else None


def save(conn, key, result, images, ttl):
    conn.execute(
        """insert into spherex.results (key, version, payload, expires_at)
           values (%s, %s, %s, now() + %s * interval '1 second')
           on conflict (key) do update set version=excluded.version,
               payload=excluded.payload, expires_at=excluded.expires_at""",
        (key, cache.VERSION, Jsonb(result), ttl),
    )
    conn.execute("delete from spherex.images where result_key=%s", (key,))
    with conn.cursor() as cur:
        cur.executemany("insert into spherex.images values (%s, %s, %s)",
                        [(key, name, data) for name, data in images])


def run(job, runner):
    # Finish work within the invocation. Postgres releases this lock on disconnect.
    with connect() as conn:
        result = read(conn, job.key)
        if result is not None:
            return 200, result
        locked, = conn.execute("select pg_try_advisory_xact_lock(hashtextextended(%s, 0))",
                               (job.key,)).fetchone()
        if not locked:
            return 202, {"status": "running", "stage": "process", "message": "Обработка"}
        result = read(conn, job.key)
        if result is not None:
            return 200, result
        with worker_slot() as workspace:
            if enabled() and workspace is None:
                return 202, {"status": "running", "stage": "queued", "message": "Ожидание обработки"}
            job.started = True
            result = runner(job)
            urls = [f[field] for f in result.get("frames", []) for field in ("image_url", "flags_url") if f.get(field)]
            if result.get("trail_url"):
                urls.append(result["trail_url"])
            directory = workspace.path / "render" if workspace is not None else cache.CACHE / job.key
            images = [(Path(url).name, (directory / Path(url).name).read_bytes()) for url in urls]
            save(conn, job.key, result, images, cache.ttl(result))
            # Commit before workspace cleanup: response images must already be durable.
            conn.commit()
        return 200, result
