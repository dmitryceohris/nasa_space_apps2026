# SNAPSHOT — Единый воспроизводимый снимок проекта SPHEREx Sky Viewer

## 1. Идентификаторы веток и коммитов (Git)

- Репозиторий: `/Users/asanalialdungarov/Desktop/nasa_space_apps2026` (`https://github.com/dmitryceohris/nasa_space_apps2026.git`)
- Рабочая ветка снимка: `asanali/spectrum-fixes`
- Базовый коммит `origin/asanali/spectrum-fixes`: `51d6941bc5b7109af914561eb430b50dc2ffbec2`
- Базовый коммит `origin/main` (включая PR #4 `fix/frontend-api`): `3be3e325eb5331dacb774cc647f48ff35828d493`
- Коммит слияния `origin/main` в `asanali/spectrum-fixes`: `a04ccb92b1e6591c82efa2a61a1a39dc7641fb9e`
- Ветка `codex/physics-five-targets`: оставлена нетронутой на `47d0d84b2376e377db3ad9382e23fe79106b2454` (в снимок НЕ входит)

## 2. Окружение и зафиксированные зависимости

- Версия Python (`.python-version`): `3.12.13` (`requires-python = "==3.12.*"`)
- Файл блокировки зависимостей: `uv.lock` (150 пакетов)
- SHA-256 файла `uv.lock`:
  `e8cc9fec06601f3f316d4f2b30c6b1ba747c2f3935ed420226e7d1dcb5ef72d7`
- Ключевые зафиксированные версии в `pyproject.toml` и `uv.lock`:
  - `astropy==8.0.2`
  - `astroquery==0.4.11`
  - `fastapi==0.162.7`
  - `httpx==0.28.1`
  - `numpy==2.5.4`
  - `pillow==13.0.0`
  - `pytest==9.2.0`
  - `reproject==0.21.0`
  - `scipy==1.18.1`
  - `uvicorn==0.54.0`

## 3. Состав снимка по зонам ответственности

- **Нурислам** (`src/Exposure/`, `src/processing/`, `api/centroid.py`): загрузка FITS SPHEREx, вычитание зодиакального света, битовые маски `FLAGS`, коаддиция по эпохам, разностные кадры, построение `trail`, апертурная фотометрия спектра и кривой блеска, детекция центроидов.
- **Дилайла и Айгерим** (`OBJECT_PARAMETERS.csv`, `OBJECT_PARAMETERS.md`, `AIGERIM_DATA_CHECK.md`, `BARNARD_MOTION.md`, `sky_motion.py`, `barnard_motion.py`, `docs/PHYSICS.md`, `frontend/index.html`): выверенные астрометрические параметры 5 целей (Gaia DR3 и Kirkpatrick et al. 2021), физическая модель собственного движения и параллакса, пользовательский веб-интерфейс на `/`.
- **Асанали** (`api/`, `docs/API.md`, `SNAPSHOT.md`, `start.sh`, `start_demo.sh`, `scripts/prepare_demo_cache.py`, `mock/`, `tests/test_api*.py`):
  - FastAPI-бэкенд (`/api/health`, `/api/targets`, `/api/sources`, `/api/resolve`, `/api/view`, `/api/export.gif`, `/api/spectrum`, `/api/lightcurve`, раздача `/` и `/dev/`);
  - строгая дедупликация параллельных запросов (`api/jobs.py`, `api/irsa.py`, `api/resolve.py`);
  - ограниченная очередь активных задач (`JOB_MAX_ACTIVE`, HTTP `429 queue_full` с заголовком `Retry-After`) и гарантированное освобождение слотов в `finally`;
  - атомарная запись всех файлов кэша (`.tmp_<uuid>_<name>` + `replace`), включая PNG, `result.json`, `exp_*.json`, `dl_*.txt`, `cut_*.fits`, `resolve.json`;
  - запрет вечного кэширования частичных сбоев (`failed > 0`), TTL для пустых ответов (`CACHE_EMPTY_TTL`) и списков экспозиций (`IRSA_EXP_TTL`);
  - экспоненциальный backoff (`IRSA_RETRIES`, `IRSA_BACKOFF`, `IRSA_BASE_URL`, `SESAME_URLS`) и возврат устаревшего кэша с `"stale": true` и предупреждением в `"warnings"` при недоступности архива (`502 upstream`).

## 4. Воспроизводимая инструкция проверки с нуля (из чистого клона)

```bash
# 1. Проверить целостность lock-файла и синхронизировать окружение
uv lock --check
uv sync --locked

# 2. Запустить полный набор офлайн-тестов (118 passed)
uv run --offline --locked pytest -q tests

# 3. Запустить сервер одной командой (автоматически создаёт .venv при отсутствии)
PREWARM=0 PORT=8000 bash start.sh
# или в демо-режиме с фиксированной эпохой CENTER_MJD=61358:
PORT=8000 bash start_demo.sh

# 4. Проверить работоспособность API
curl -i http://127.0.0.1:8000/api/health
curl -s http://127.0.0.1:8000/api/targets
curl -s "http://127.0.0.1:8000/api/resolve?target=barnard"
```
