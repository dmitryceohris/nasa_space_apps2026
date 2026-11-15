"""Export the five challenge targets through the real /api/view and GIF routes.

Replays the saved archive metadata and exact local FITS, with all network access
disabled. Processing, calibration, epoch grouping, masks, render and HTTP routes
are the application code. Writes a public gallery at /challenge-demo/.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import logging
import os
from pathlib import Path
import sys
import tempfile
import time
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]
os.environ["PREWARM"] = "0"

from astropy.io import fits
from astropy.time import Time
from astropy.utils import iers
from fastapi.testclient import TestClient
from PIL import Image, ImageDraw, ImageFont

from api import cache, irsa, main, resolve
from api.limits import RequestLimiter
from processing.epochs import select_per_epoch

SNAPSHOTS = {
    "barnard": "08e6182e8643a70c",
    "epsindib": "7906fc5eb012df88",
    "wise0855": "19aace5c498b685a",
    "rwaur": "66f189f43ec22f35",
    "rdra": "3face6ff239f7a34",
}
NAMES = {"barnard": "Звезда Барнарда", "epsindib": "Epsilon Indi B",
         "wise0855": "WISE 0855", "rwaur": "RW Aurigae", "rdra": "R Draconis"}
AUDIT = ROOT / "data/validation/science-remediation/measurements.json"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dump(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n")


def counted(value, forms):
    form = forms[2] if 11 <= value % 100 <= 14 else forms[0] if value % 10 == 1 else forms[1] if value % 10 in (2, 3, 4) else forms[2]
    return f"{value} {form}"


def no_network(*args, **kwargs):
    raise AssertionError("This export uses saved real inputs only; network disabled")


def font(size, bold=False):
    folder = Path("/System/Library/Fonts/Supplemental")
    paths = [folder / ("Arial Bold.ttf" if bold else "Arial.ttf"),
             Path("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf")]
    for path in paths:
        if path.exists():
            return ImageFont.truetype(str(path), size)
    raise FileNotFoundError("Install a Cyrillic font or edit font() with its path")


def annotate(response, raw_frames, out, name):
    """Add metadata only. Sky pixels retain nearest-neighbour sampling and order."""
    images = []
    for i, (frame, path) in enumerate(zip(response["frames"], raw_frames)):
        sky = Image.open(path).convert("RGB")
        canvas = Image.new("RGB", (720, 940), "#101722")
        d = ImageDraw.Draw(canvas)
        d.text((40, 24), name, font=font(32, True), fill="white")
        d.text((40, 70), f"{frame['date']} UTC    Кадр {i+1}/{len(raw_frames)}", font=font(24), fill="#dce5f2")
        wave = f"{frame['wave_um']:.3f}" if frame["wave_um"] is not None else "неизвестна"
        d.text((40, 108), f"D{frame['det']}    λ = {wave} мкм    Экспозиций: {frame['n_exp']}",
               font=font(22), fill="#c4d2e5")
        d.text((40, 145), "Север сверху, восток слева", font=font(19), fill="#94a6bf")
        # The 640-pixel square is a display enlargement, not extra resolution.
        canvas.paste(sky.resize((640, 640), Image.Resampling.NEAREST), (40, 180))
        pred = frame.get("predicted_px")
        if pred is not None:
            x = 40 + (pred[0] + .5) * 640 / sky.width
            y = 180 + (pred[1] + .5) * 640 / sky.height
            d.ellipse((x-19, y-19, x+19, y+19), outline="#55efc4", width=3)
        label = "Кольцо: прогноз по каталогу" if pred is not None else "Поле по заданным RA / Dec"
        d.text((40, 831), label, font=font(19), fill="#55efc4")
        d.text((40, 861), "NASA SPHEREx / IRSA    Фиолетовый: нет данных", font=font(18), fill="#a7b7cd")
        d.text((40, 901), "Разные полосы. Движение и переменность не подтверждены.",
               font=font(18), fill="#ffe1a8")
        images.append(canvas)
    images[0].save(out, save_all=True, append_images=images[1:], duration=900, loop=0,
                   disposal=2, optimize=False)
    with Image.open(out) as gif:
        assert gif.n_frames == len(raw_frames), (out, gif.n_frames, len(raw_frames))
        assert gif.size == (720, 940) and gif.info["loop"] == 0


def saved_input(tid):
    snapshot = ROOT / "cache" / SNAPSHOTS[tid] / "result.json"
    old = json.loads(snapshot.read_text())
    q = old["query"]
    metadata_path = ROOT / "cache" / f"exp_{q['ra']:.3f}_{q['dec']:.3f}_d{q['det']}.json"
    saved = json.loads(metadata_path.read_text())
    rows = saved["rows"]
    selected = [rows[i] for group in select_per_epoch([r["mjd"] for r in rows], 6) for i in group]
    lookup, manifest = {}, []
    for row in selected:
        uri = row["uri"]
        variants = [uri, uri.removeprefix(irsa.IRSA + "/")]
        candidates = [ROOT / "cache" / ("cut_" + hashlib.md5(
            f"{v}{q['ra']:.4f}{q['dec']:.4f}{q['fov']/60*1.5:.3f}".encode()).hexdigest() + ".fits")
            for v in variants]
        path = next((p for p in candidates if p.exists()), None)
        if path is None:
            raise FileNotFoundError(f"Exact saved FITS missing for {tid}: {uri}")
        with fits.open(path) as hdul:
            hdr = hdul["IMAGE"].header
            assert int(hdr["DETECTOR"]) == q["det"]
            assert str(hdr["OBSID"]) == irsa.obs_id(uri)
            manifest.append({"obs_id": str(hdr["OBSID"]), "detector": int(hdr["DETECTOR"]),
                             "mjd_utc": float(hdr["MJD-AVG"]), "date_utc": hdr["DATE-AVG"],
                             "product_url": uri, "path": str(path.relative_to(ROOT)), "sha256": sha(path)})
        lookup[uri] = path
    return old, rows, lookup, manifest, {"metadata_path": str(metadata_path.relative_to(ROOT)),
                                       "metadata_sha256": sha(metadata_path),
                                       "snapshot_path": str(snapshot.relative_to(ROOT)), "snapshot_sha256": sha(snapshot)}


def poll(client, query):
    started = time.monotonic()
    while True:
        response = client.get("/api/view", params=query)
        if response.status_code != 202:
            assert response.status_code == 200, response.text
            result = response.json()
            assert result["status"] == "ready" and result["frames"], result
            return result
        if time.monotonic() - started > 120:
            raise TimeoutError(query)
        time.sleep(.025)


def export_one(client, tid, query, rows, lookup, inputs, source, dest):
    dest.mkdir(parents=True, exist_ok=True)
    replay_warning = "Сохранённая выборка IRSA, повторно обработанная локально; новый поиск архива не выполнялся"
    metadata = irsa.ExposureRows(rows, source="saved_sia2_replay", warnings=[replay_warning])
    with patch.object(irsa, "exposures", lambda *args: metadata), \
         patch.object(irsa, "cutout", lambda uri, *args: lookup[uri]):
        result = poll(client, query)
        raw = client.get("/api/export.gif", params=query)
        assert raw.status_code == 200 and raw.headers["content-type"].startswith("image/gif")
        (dest / "api.gif").write_bytes(raw.content)
        dump(dest / "api_response.json", result)
        frame_paths = []
        for i, frame in enumerate(result["frames"]):
            image = client.get(frame["image_url"])
            assert image.status_code == 200 and image.content.startswith(b"\x89PNG")
            path = dest / f"frame_{i:02}.png"
            path.write_bytes(image.content)
            frame_paths.append(path)
            flags = client.get(frame["flags_url"])
            assert flags.status_code == 200
            (dest / f"frame_{i:02}_flags.png").write_bytes(flags.content)
    assert sum(f["n_exp"] for f in result["frames"]) == result["exposure_counts"]["used"]
    assert result["exposure_counts"]["failed"] == 0
    assert len({(x["obs_id"], x["detector"]) for x in inputs}) == len(inputs)
    assert {o for f in result["frames"] for o in f["obs_ids"]} == {x["obs_id"] for x in inputs}
    annotate(result, frame_paths, dest / "blink.gif", NAMES[tid])
    with Image.open(dest / "api.gif") as native:
        assert native.n_frames == len(result["frames"])
    dump(dest / "inputs.json", {"source": source, "inputs": inputs})
    frames = result["frames"]
    waves = [f["wave_um"] for f in frames if f["wave_um"] is not None]
    mjds = [a["mjd_utc"] for a in inputs]
    stats = {"id": tid, "name": NAMES[tid], "detector": query["det"], "fov_arcmin": query["fov"],
             "mode": query["mode"], "frames": len(frames), "exposures": result["exposure_counts"]["used"],
             "exposure_counts": result["exposure_counts"],
             "first_exposure_utc": Time(min(mjds), format="mjd").isot,
             "last_exposure_utc": Time(max(mjds), format="mjd").isot,
             "exposure_baseline_days": round(max(mjds)-min(mjds), 6),
             "first_frame_date": frames[0]["date"], "last_frame_date": frames[-1]["date"],
             "frame_wavelength_range_um": [min(waves), max(waves)] if waves else None,
             "min_coverage_fraction": min(f["quality"]["coverage_fraction"] for f in frames),
             "max_coverage_fraction": max(f["quality"]["coverage_fraction"] for f in frames),
             "candidate_display_markers": sum(f["measured_px"] is not None for f in frames),
             "validated_astrometry_frames": sum(f["validation"]["status"] == "validated" for f in frames),
             "measured_motion_arcsec_yr": result.get("motion", {}).get("pm_measured_arcsec_yr") if result.get("motion") else None,
             "warnings": result["warnings"], "gif": f"{dest.name}/blink.gif",
             "api_gif": f"{dest.name}/api.gif", "api_response": f"{dest.name}/api_response.json"}
    dump(dest / "statistics.json", stats)
    print(json.dumps({k: stats[k] for k in ("id", "frames", "exposures", "first_frame_date", "last_frame_date")}), flush=True)
    return stats


def science_statistics():
    data = json.loads(AUDIT.read_text())
    out = []
    for tid in SNAPSHOTS:
        rows = [m for m in data["measurements"] if m["target"] == tid]
        products = defaultdict(set)
        for r in rows:
            products[r["product_url"]].add(r["status"])
        # Do not silently count a product with contradictory statuses twice.
        assert all(len(v) == 1 for v in products.values())
        counts = Counter(next(iter(v)) for v in products.values())
        out.append({"id": tid, "name": NAMES[tid], "unique_products": len(products),
                    "accepted": counts["validated"], "rejected": counts["rejected"],
                    "unavailable": counts["unavailable"], "attempts": len(rows),
                    "validated_motion": data["summary"][tid]["validated_motion"],
                    "validated_lightcurve": data["summary"][tid]["validated_lightcurve"]})
    return {"source": str(AUDIT.relative_to(ROOT)), "source_sha256": sha(AUDIT),
            "checked_at_utc": data["checked_at_utc"], "scope": "Separate multi-detector science audit; not the GIF exposure subset",
            "deduplication": "Full physical product_url, including release and detector; repeated cutouts counted once",
            "targets": out}


def write_reports(out, targets, objective, science):
    total_frames = sum(t["frames"] for t in targets)
    total_exps = sum(t["exposures"] for t in targets)
    report = {"created_at_utc": Time.now().isot, "source": "Application /api/view and /api/export.gif with saved real archive inputs",
              "network_requests": 0, "cache_version": cache.VERSION, "center_mjd": 61358,
              "credits": "NASA SPHEREx via NASA/IPAC IRSA", "total_frames": total_frames,
              "total_exposures": total_exps, "targets": targets, "wise_objective_coordinate_test": objective,
              "science_audit": science,
              "source_sha256": {str(p.relative_to(ROOT)): sha(p) for folder in ("api", "src") for p in (ROOT/folder).rglob("*.py")},
              "limitations": ["Saved bounded archive sample, not a fresh archive-wide census",
                              "GIF ring is a catalogue prediction, not an independently measured position",
                              "Coadds can combine different narrow spectral channels; brightness changes do not establish variability",
                              "Native sky resolution is 6.15 arcsec/pixel; GIFs are enlarged with nearest-neighbour sampling",
                              "Validated photometry does not validate astrometry, motion, or a light curve"]}
    dump(out / "statistics.json", report)
    fields = ["id", "name", "frames", "exposures", "detector", "fov_arcmin", "first_exposure_utc",
              "last_exposure_utc", "exposure_baseline_days", "first_frame_date", "last_frame_date",
              "min_coverage_fraction", "max_coverage_fraction", "validated_astrometry_frames"]
    with (out / "statistics.csv").open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(targets)
    table = "\n".join(f"| {t['name']} | {t['frames']} | {t['exposures']} | D{t['detector']} | {t['first_exposure_utc'][:10]} — {t['last_exposure_utc'][:10]} | {t['exposure_baseline_days']:.1f} |" for t in targets)
    science_table = "\n".join(f"| {t['name']} | {t['unique_products']} | {t['accepted']} | {t['rejected']} | {t['unavailable']} |" for t in science["targets"])
    md = f"""# GIF и статистика пяти целей SPHEREx

Повторная обработка сохранённых реальных FITS текущим кодом приложения, 15 ноября 2026 года.
Новый поиск IRSA не выполнялся. Все GIF получены через `/api/view` и `/api/export.gif`.

## Просмотр приложения

| Цель | Кадров GIF | Экспозиций | Детектор | Даты исходных экспозиций, UTC | База, сут |
|---|---:|---:|---|---|---:|
{table}
| **Всего** | **{total_frames}** | **{total_exps}** | | | |

Число кадров означает число эпох после группировки. Один кадр может объединять несколько экспозиций.
Дата в GIF — средняя дата кадра, границы таблицы — даты крайних исходных экспозиций.
Кольцо показывает прогноз каталожной модели. Поле каждого GIF фиксировано между эпохами.
Одинаковая шкала яркости действует на всю последовательность, без интерполяции между датами.
Фиолетовые пиксели обозначают отсутствие пригодных данных. Увеличение не добавляет разрешения.
Разные спектральные каналы и маски могут менять внешний вид. Подтверждённых астрометрических кадров: 0.

## Отдельная научная проверка

Эта проверка включает другие экспозиции и все детекторы. Её объём не равен объёму GIF.
Источник: `{science['source']}`, проверка `{science['checked_at_utc']}` UTC.
Повторные вырезки одного продукта исключены по полному URL продукта.

| Цель | Уникальных продуктов | Принятая фотометрия | Отклонено | Недоступно |
|---|---:|---:|---:|---:|
{science_table}
| **Всего** | **214** | **6** | **188** | **20** |

Принятая фотометрия WISE 0855 выполнена в каталожно заданной позиции.
Собственное движение и переменность ни одной из пяти целей в этом наборе не подтверждены.
Числа относятся к ограниченным сохранённым выборкам, а не ко всему архиву SPHEREx.

## Проверка точного ввода координат WISE 0855

RA `08h55m10.831s`, Dec `-7d14m42.53s` переданы в `/api/view` без имени цели.
Дополнительная последовательность `wise0855-coordinates/blink.gif`: {objective['frames']} кадра,
{objective['exposures']} экспозиции, D6, поле 6′. Это фиксированное поле без модели движения.
Она не включена в итог {total_frames} кадров и {total_exps} экспозиций пяти основных GIF.

## Файлы

- `index.html`: общая галерея, доступная по `/challenge-demo/` при запуске приложения.
- `<target>/blink.gif`: версия с названием, датой, длиной волны и каталожным прогнозом.
- `<target>/api.gif`: исходный GIF-экспорт приложения без подписей и увеличения.
- `<target>/api_response.json`: полный ответ приложения. Его `/cache/` URL относятся к изолированному прогону.
- `<target>/frame_*.png`: локальные копии кадров и масок для просмотра без временного кэша.
- `<target>/inputs.json`: исходные продукты, контрольные суммы FITS и сохранённых метаданных.
- `statistics.json`, `statistics.csv`: статистика приложения и отдельная научная проверка.

Воспроизведение из корня проекта: `.venv/bin/python scripts/export_challenge_gifs.py`.
Нужен исходный локальный FITS-кэш, указанный в манифестах. Обращения к сети запрещены самим экспортом.
"""
    (out / "STATISTICS.md").write_text(md)
    cards = "\n".join(f'''<article><h2>{t['name']}</h2><img src="{t['gif']}" alt="Последовательность SPHEREx: {t['name']}"><p>{counted(t['frames'], ('кадр','кадра','кадров'))} из {counted(t['exposures'], ('экспозиции','экспозиций','экспозиций'))}, D{t['detector']}<br>{t['first_exposure_utc'][:10]} — {t['last_exposure_utc'][:10]} UTC</p><a href="{t['gif']}" download>Скачать GIF</a> · <a href="{t['api_gif']}" download>Исходный GIF приложения</a></article>''' for t in targets)
    rows = "\n".join(f"<tr><td>{t['name']}</td><td>{t['unique_products']}</td><td>{t['accepted']}</td><td>{t['rejected']}</td><td>{t['unavailable']}</td></tr>" for t in science["targets"])
    html = f'''<!doctype html><html lang="ru"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><link rel="icon" href="data:,"><title>Пять целей SPHEREx</title>
<style>body{{margin:0;background:#101722;color:#e6edf8;font:17px/1.6 system-ui,sans-serif}}main{{max-width:1280px;margin:auto;padding:40px 24px}}h1{{font-size:38px;line-height:1.2}}h2{{font-size:24px}}a{{color:#55efc4}}.gallery{{display:grid;grid-template-columns:repeat(auto-fit,minmax(300px,1fr));gap:36px}}article img{{width:100%;height:auto}}p{{max-width:1000px;color:#bccbdf}}table{{border-collapse:collapse;width:100%;margin:24px 0}}td,th{{padding:12px;text-align:left;border-bottom:1px solid #344255}}.scroll{{overflow:auto}}.note{{color:#ffe1a8}}footer{{margin-top:50px;color:#94a6bf}}</style>
<main><h1>Пять целей SPHEREx во времени</h1><p>{counted(total_frames, ('кадр','кадра','кадров'))} из {counted(total_exps, ('экспозиции','экспозиций','экспозиций'))}. Повторная обработка сохранённых реальных FITS кодом приложения от 15 ноября 2026 года.</p>
<p class="note">Кольцо показывает каталожный прогноз положения. Движение и переменность по этим GIF не подтверждены. Спектральные полосы меняются между кадрами. Фиолетовый обозначает отсутствие пригодных данных.</p>
<p><a href="statistics.csv" download>Статистика CSV</a> · <a href="statistics.json" download>Полные данные JSON</a> · <a href="STATISTICS.md">Научная записка</a></p>
<section class="gallery">{cards}</section><h2>Отдельная научная проверка</h2><p>Другая, более широкая выборка по всем детекторам. Повторные вырезки исключены. Принятая фотометрия относится к измерению потока в каталожной позиции и не подтверждает движение.</p>
<div class="scroll"><table><thead><tr><th>Цель</th><th>Продуктов</th><th>Принято</th><th>Отклонено</th><th>Недоступно</th></tr></thead><tbody>{rows}<tr><th>Всего</th><th>214</th><th>6</th><th>188</th><th>20</th></tr></tbody></table></div>
<h2>Проверка ввода RA и Dec</h2><p>Для WISE 0855 отдельно проверены координаты из задания: RA 08h55m10.831s, Dec −7d14m42.53s. Фиксированное поле 6′, D6, {objective['frames']} кадра. <a href="wise0855-coordinates/blink.gif">Открыть GIF координатного запроса</a>.</p>
<footer>NASA SPHEREx via NASA/IPAC IRSA. Сохранённые реальные FITS, новый поиск архива не выполнялся.</footer></main></html>'''
    (out / "index.html").write_text(html)
    return report


def run(out):
    out.mkdir(parents=True, exist_ok=True)
    logging.getLogger("httpx").setLevel(logging.WARNING)
    iers.conf.auto_download = False
    targets = []
    with tempfile.TemporaryDirectory(prefix="spherex-gif-export-") as temporary:
        private = Path(temporary)
        mount = next(r.app for r in main.app.routes if getattr(r, "name", "") == "cache")
        with patch.object(cache, "CACHE", private), \
             patch.object(mount, "all_directories", [str(private)]), \
             patch.object(main, "request_limiter", RequestLimiter(100000)), \
             patch.dict(os.environ, {"CENTER_MJD": "61358", "API_RATE_LIMIT_PER_MINUTE": "100000"}), \
             patch.object(irsa, "_send", no_network), patch("httpx.request", no_network), \
             patch("httpx.get", no_network):
            with TestClient(main.app) as client:
                for tid in SNAPSHOTS:
                    old, rows, lookup, inputs, source = saved_input(tid)
                    q = {"target": tid, "fov": old["query"]["fov"], "det": old["query"]["det"],
                         "mode": "coadd", "max_per_epoch": 6}
                    targets.append(export_one(client, tid, q, rows, lookup, inputs, source, out / tid))
                evidence_path = ROOT / "data/validation/science-acceptance/evidence.json"
                evidence = json.loads(evidence_path.read_text())["wise_exact_coordinate_frames"]
                needed = {r["cutout_sha256"] for r in evidence}
                paths = {sha(p): p for p in (ROOT/"cache").glob("cut_*.fits")}
                assert needed <= paths.keys()
                rows = [{"uri": r["product_url"], "mjd": r["fits_header"]["MJD-AVG"]} for r in evidence]
                lookup = {r["product_url"]: paths[r["cutout_sha256"]] for r in evidence}
                inputs = [{"obs_id": r["obs_id"], "detector": r["detector"],
                           "mjd_utc": r["fits_header"]["MJD-AVG"], "date_utc": r["date_utc"],
                           "product_url": r["product_url"], "sha256": r["cutout_sha256"],
                           "path": str(paths[r["cutout_sha256"]].relative_to(ROOT))} for r in evidence]
                q = {"ra": "08h55m10.831s", "dec": "-7d14m42.53s", "fov": 6, "det": 6,
                     "mode": "single", "max_per_epoch": 1}
                objective = export_one(client, "wise0855", q, rows, lookup, inputs,
                                       {"evidence_path": str(evidence_path.relative_to(ROOT)), "sha256": sha(evidence_path)},
                                       out / "wise0855-coordinates")
                assert client.get("/api/health").status_code == 200
    report = write_reports(out, targets, objective, science_statistics())
    # Verify the delivered static paths using the unmodified application mount.
    with TestClient(main.app) as client:
        assert client.get("/challenge-demo/").status_code == 200
        assert client.get("/challenge-demo/statistics.json").json()["total_frames"] == report["total_frames"]
        for t in targets:
            assert client.get("/challenge-demo/" + t["gif"]).content == (out/t["gif"]).read_bytes()
    print(json.dumps({"output": str(out), "total_frames": report["total_frames"],
                      "total_exposures": report["total_exposures"], "network_requests": 0}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "frontend/challenge-demo")
    args = parser.parse_args()
    run(args.output.resolve())
