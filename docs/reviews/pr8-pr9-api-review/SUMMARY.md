# Итог этапа: API-ревью PR #8 и PR #9 (2026-11-14)

Исторический отчёт о проверке указанных ниже коммитов. Последующие решения,
исправленная альтернатива закрытия #8 и новые проверки описаны в
[отчёте интеграции](../asanali-integration-2026-11-14.md).

Ревьюер API-части: Асанали. Репозиторий: dmitryceohris/nasa_space_apps2026. Тексты опубликованных комментариев лежат рядом.

## Головы веток на момент публикации

| Ветка | Голова | Состояние |
|---|---|---|
| origin/main | afa35ef | – |
| PR #8 `codex/processing-acceptance` | df04803 (было c435817 до 23:04 +0500) | OPEN |
| PR #9 `codex/science-final-acceptance` | 03c22bc | OPEN |
| PR #7 `asanali/reliability-main` | e93d6a3 | CLOSED, ветка сохранена |

## Публикации

| PR | Вид | Время (UTC) | Ссылка |
|---|---|---|---|
| #7 | комментарий + закрытие (без удаления ветки) | 2026-11-14 18:12:44 | https://github.com/dmitryceohris/nasa_space_apps2026/pull/7#issuecomment-6882022587 |
| #9 | review, state COMMENTED (не APPROVED), «согласовано с ограничениями» | 2026-11-14 18:18:22 | https://github.com/dmitryceohris/nasa_space_apps2026/pull/9#pullrequestreview-5906296530 |
| #8 | комментарий | 2026-11-14 18:18:45 | https://github.com/dmitryceohris/nasa_space_apps2026/pull/8#issuecomment-6882101154 |

Тексты: `pr7_close_comment.md`, `pr9_review_comment_v2.md`, `pr8_comment_v2.md` (в этой папке). Файлов v1 в этой ветке нет.

## Тесты (свежие клоны, `uv run --offline --locked pytest -q`)

| Ветка | Результат | Прогонов |
|---|---|---|
| main afa35ef | 225 passed, 1 skipped | 1 |
| #8 c435817 | 242 passed, 1 skipped | 3 |
| #8 df04803 | 248 passed, 1 skipped | 3 |
| #9 03c22bc | 259 passed, 1 skipped | 3 |

Skip во всех случаях: `tests/test_exposure.py:68` (нужен `SPHEREX_LIVE=1`).
`tests/test_api_acceptance.py` + `tests/test_api_public_limits.py`: 5 + 5 passed на #9 и на #8 df04803.

Docker (#9): образ `spherex-api:pr9`, контейнер `spherex-pr9` (порт 8050), том `spherex-pr9-cache` оставлены. `/api/health` 200, `/api/resolve` без target 400 `bad_param`, `/api/no-such` 404 `not_found`. Недоступный архив: 202, затем 200 `status: error`, тот же `job_id` при повторе, ошибка не кэшируется, 3 попытки к IRSA с паузами 5 и 10 с.

## Не проверено

- Живой `/api/view` на spherex-pr9 (wise0855 и фиксированные RA/Dec до `ready`, `target_association` и `offset_from_prediction_px` на реальных кадрах, `processing.registration`). Причина: IRSA TAP 502 и таймауты с 17:42:27 по 17:52:25 UTC; рабочий TAP-запрос из контейнера завершался таймаутом также в 17:56:43, 17:58:52 и 18:20:02 UTC.
- Прогрев демо-кэша после подъёма `cache.VERSION`: прежний кэш не отдаётся, перед показом пересобрать `scripts/prepare_demo_cache.py`; время не замерено.

## Рекомендуемый порядок слияния

1. Слить #9 первым как более полный.
2. Автору #8 перебазировать df04803 на main после слияния #9 и оставить только то, чего нет в #9 (`prediction_outside_field`, отказ от повторов на 4xx в `api/irsa.py`, удаление `Job.result`/`ready`, правки вне API-зоны), с `cache.VERSION = 11`. Альтернатива: перенести в #9 все нужные уникальные изменения #8 с их тестами. Двух API-исправлений недостаточно: нужно также учесть ассоциацию/assisted-измерения и переход RA 0/360 в `spxhack`, все отрицательные доли в centroid, отсутствие VARIANCE в loader и frontend-запрос спектра. Закрывать #8 только после проверки каждого уникального изменения; сохранять защиты SOURCE, цвет отсутствующих данных и научные ограничения из #9.

Пробное слияние симметрично: оба порядка дают конфликты в одних и тех же 10 файлах (21 блок); в API-зоне `api/cache.py` (1), `api/jobs.py` (4), `docs/API.md` (4).

## Открытые вопросы к авторам

| Вопрос | Кому | Где |
|---|---|---|
| Документация `spectral_sampling`: не описан `n_samples` и то, что без отсчётов все поля `null` (включая `different_channels`); у кривой блеска другой состав (нет `n_samples`) и нет округления | #9 | `docs/API.md` ~790, `src/processing/frames.py:196-212`, `src/processing/pipeline.py:182` |
| Три постоянных warnings (около 370 символов) в каждом ответе с кадрами; предложено перенести в `processing.registration.limitations` | #9 (решает автор) | `api/jobs.py:278-288`, `frontend/index.html:175`, `:239-240` |
| `prediction_outside_field` не обрабатывается во frontend: при `status ok` подпись «centroid candidate…», хотя маркер скрыт | #8 (после перебазирования) | `frontend/index.html`, `measuredText` |
| `cache.VERSION` до 11 при перебазировании #8; в #8 `docs/API.md` противоречит себе (стр. ~766 «8», «Дополнение» «9»); df04803 менял `result.json` без подъёма версии | #8 | `api/cache.py:17`, `docs/API.md` |
| Тест, что 408/429/5xx по-прежнему повторяются (сейчас тест только на отсутствие повторов на 404) | #8 | `tests/test_processing_acceptance.py`, `api/irsa.py` |

Научная корректность порогов (`ASSOC_PX` и др.), масок и фотометрии API-ревью не оценивалась.
