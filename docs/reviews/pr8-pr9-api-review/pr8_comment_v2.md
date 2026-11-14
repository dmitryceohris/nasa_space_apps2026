Спасибо за работу над B1 и за быстрые правки по ревью. Ниже заметка по API-части (ветка `df04803`) и предложение по порядку слияния.

### Что проверено на `df04803`

- Свежий клон, `uv run --offline --locked pytest -q`, три прогона подряд: 248 passed, 1 skipped в каждом (skip: `tests/test_exposure.py:68`, нужен `SPHEREX_LIVE=1`). `tests/test_api_acceptance.py` и `tests/test_api_public_limits.py`: 5 + 5 passed.
- Фиксированные RA/Dec через `/api/view` (офлайн, синтетика): 202, затем 200 `ready`; у каждого кадра `prediction_status = no_model`, `target_association = no_prediction`, `offset_from_prediction_px = null`.
- `prediction_outside_field` описан в `docs/API.md` и покрыт тестом `test_prediction_outside_the_field_gives_no_measured_position`.
- Удаление `Job.result` и статуса `ready` безопасно: на `main` они нигде не присваивались (`j.status` получает только `"error"`), других ссылок в `api/`, `tests/` и `frontend/` нет. Поведение 202 -> 200 не изменилось.
- `api/irsa.py`: 400, 403 и 404 дают одну попытку и 502 `upstream` «IRSA отклонил запрос», а 408, 429, 500, 502 и 503 по-прежнему дают 3 попытки (проверено подменой `httpx.request`). В репозитории тест есть только на 404.

### Связь с #9

PR #9 перенёс прежнюю версию этого PR (`246603f` -> `743c6e7` с правками) и добавляет `processing.registration` и `cache.VERSION = 10`. Новый коммит `df04803` в #9 не входит. Только здесь есть:

- `target_association = "prediction_outside_field"` и `measured_px = null` в этом случае (`api/jobs.py`, `docs/API.md`);
- отказ от повторов на 4xx в `api/irsa.py` и тест на него;
- удаление неиспользуемых `Job.result` и статуса `ready`;
- вне API-зоны: правки `src/processing/` (centroid, frames, photometry, render), `src/Exposure/loader.py`, `src/spxhack/measure.py`, `frontend/index.html` (запрос спектра без `per_epoch`), README, тесты и отчёт `docs/reviews/processing-acceptance.md`.

Пробное слияние в свежем клоне симметрично. И «#9, затем #8», и «#8, затем #9» дают конфликты в одних и тех же 10 файлах (21 блок):
`api/cache.py`, `api/jobs.py`, `docs/API.md`, `docs/reviews/processing-acceptance.md`, `src/processing/centroid.py`, `src/processing/frames.py`, `src/processing/photometry.py`, `src/processing/pipeline.py`, `src/processing/render.py`, `tests/test_processing_acceptance.py`.

### Предложение

Слить #9 первым как более полный. Затем перебазировать эту ветку на `main` и оставить в PR только то, чего нет в #9 (список выше). Так как меняется содержимое `result.json` (новое значение `target_association`), поднять `cache.VERSION` до 11 и обновить номер версии в `docs/API.md`.

### Замечания для перебазированной версии

1. `frontend/index.html` не читает `target_association` (в этой ветке нигде; в #9 обрабатывается только `far_from_prediction`). При `prediction_outside_field` с `measurement.status = ok` подпись кадра будет «centroid candidate: good ±… px», хотя маркер измерения скрыт и показаны кандидаты; рядом в `warnings` будет «Эпоха N: прогноз вне поля; маркер скрыт». Нужна одна ветка в `measuredText` по аналогии с `far_from_prediction` из #9.
2. `docs/API.md` здесь противоречит сам себе: около строки 766 «Версия кэша 8 инвалидирует прежние результаты», а в разделе «Дополнение» «Версия runtime-кэша — **9**». При перебазировании оба места привести к новой версии.
3. Для сведения: `df04803` меняет содержимое `result.json` относительно `c435817`, но `cache.VERSION` в ветке остался 9. Если где-то уже есть кэш, собранный с `c435817`, он будет отдаваться со старой семантикой. При подъёме до 11 вопрос снимается.
4. Желателен тест, что 408/429/5xx по-прежнему повторяются: сейчас тест закрепляет только отсутствие повторов на 404.
