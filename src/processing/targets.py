"""Тестовые цели челленджа и пересчёт их положения на дату снимка.

Параметры берутся из OBJECT_PARAMETERS.csv (Gaia DR3 / Kirkpatrick+2021), который собрала команда.
Пересчёт — astropy `apply_space_motion` (барицентрический ICRS, без годичного параллакса:
для маркера на картинке с пикселем 6.15″ этого достаточно; параллакс Барнарда 0.55″ < 0.1 px).
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.coordinates import Distance, SkyCoord
from astropy.time import Time

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CATALOG = REPO_ROOT / "OBJECT_PARAMETERS.csv"

# id -> (имя в CSV, детектор по умолчанию, поле зрения в угл. минутах, подсказка)
PRESETS: dict[str, tuple[str, int, float, str]] = {
    "barnard": ("Barnard's Star", 3, 3.0, "самое быстрое собственное движение, ~10.4″/год"),
    "epsindib": ("Epsilon Indi B", 3, 4.0, "пара коричневых карликов: ищите метан в спектре"),
    "wise0855": ("WISE 0855-0714", 6, 3.0, "самый холодный коричневый карлик: виден только на 4–5 мкм"),
    "rwaur": ("RW Aur", 2, 3.0, "молодая переменная звезда T Tauri"),
    "rdra": ("R Dra", 2, 3.0, "мирида, период ~245 дней"),
}


@dataclass
class Target:
    id: str
    name: str
    ra: float              # град, на эпоху ref_epoch
    dec: float
    pmra: float            # mas/год (μα·cosδ)
    pmdec: float           # mas/год
    parallax: float        # mas (0 = неизвестен)
    rv: float              # км/с (0 = неизвестна)
    ref_epoch: Time
    det: int = 3
    fov: float = 3.0
    note: str = ""
    source: str = ""
    parallax_error_mas: float | None = None

    def coord_at(self, time: Time | float) -> SkyCoord:
        """Положение на момент `time` (Time или MJD)."""
        t = time if isinstance(time, Time) else Time(float(time), format="mjd", scale="utc")
        if self.pmra == 0 and self.pmdec == 0:
            return SkyCoord(self.ra * u.deg, self.dec * u.deg)
        plx = self.parallax if self.parallax > 0 else 1e-3  # далёкий объект, если параллакс неизвестен
        c = SkyCoord(ra=self.ra * u.deg, dec=self.dec * u.deg,
                     pm_ra_cosdec=self.pmra * u.mas / u.yr, pm_dec=self.pmdec * u.mas / u.yr,
                     distance=Distance(parallax=plx * u.mas),
                     radial_velocity=self.rv * u.km / u.s, obstime=self.ref_epoch)
        moved = c.apply_space_motion(new_obstime=t)
        return SkyCoord(moved.ra, moved.dec)

    def at(self, time: Time | float) -> SkyCoord:
        """Короткое имя для coord_at."""
        return self.coord_at(time)

    def to_json(self) -> dict:
        return {"id": self.id, "name": self.name, "ra": self.ra, "dec": self.dec,
                "pmra_masyr": self.pmra, "pmdec_masyr": self.pmdec,
                "parallax_mas": self.parallax, "ref_epoch": float(self.ref_epoch.jyear),
                "det": self.det, "fov": self.fov, "note": self.note, "source": self.source}


def _f(row: dict, key: str) -> float:
    v = (row.get(key) or "").strip()
    return float(v) if v else 0.0


def _epoch(row: dict) -> Time:
    value, fmt = row["ref_epoch"].strip(), (row.get("ref_epoch_format") or "jyear").strip()
    scale = (row.get("ref_epoch_scale") or "").strip().lower() or "utc"
    return Time(float(value), format=fmt, scale=scale)


def load_targets(path: str | Path = DEFAULT_CATALOG) -> dict[str, Target]:
    with open(path, newline="", encoding="utf-8") as f:
        rows = {r["object"].strip(): r for r in csv.DictReader(f)}
    out: dict[str, Target] = {}
    for tid, (name, det, fov, note) in PRESETS.items():
        r = rows.get(name)
        if r is None:
            continue
        out[tid] = Target(
            id=tid, name=name, ra=_f(r, "ra_deg"), dec=_f(r, "dec_deg"),
            pmra=_f(r, "pmra_mas_yr"), pmdec=_f(r, "pmdec_mas_yr"),
            parallax=_f(r, "parallax_mas"), rv=_f(r, "radial_velocity_km_s"),
            ref_epoch=_epoch(r), det=det, fov=fov, note=note,
            source=(r.get("catalog") or "").strip(),
            parallax_error_mas=_f(r,"parallax_error_mas") if r.get("parallax_error_mas") else None)
    return out


def find_preset(query: str, targets: dict[str, Target]) -> Target | None:
    """Ищет пресет по id или имени без учёта регистра/пробелов/апострофов."""
    norm = lambda s: "".join(ch for ch in s.lower() if ch.isalnum())  # noqa: E731
    q = norm(query)
    for t in targets.values():
        if q in (norm(t.id), norm(t.name)):
            return t
    aliases = {"barnardsstar": "barnard", "epsind": "epsindib", "epsilonindib": "epsindib",
               "wise08550714": "wise0855", "wisej08551083071442.5": "wise0855"}
    tid = aliases.get(q)
    return targets.get(tid) if tid else None


def resolve_simbad(name: str) -> Target | None:
    """Любое имя через SIMBAD (сеть). Возвращает Target с собственным движением на J2000."""
    from astroquery.simbad import Simbad

    s = Simbad()
    s.add_votable_fields("pmra", "pmdec", "plx_value")
    tbl = s.query_object(name)
    if tbl is None or len(tbl) == 0:
        return None
    row = tbl[0]
    cols = {c.lower(): c for c in tbl.colnames}

    def g(*keys):
        for k in keys:
            if k in cols:
                v = row[cols[k]]
                try:
                    return float(v) if not np.ma.is_masked(v) else 0.0
                except (TypeError, ValueError):
                    return v
        return 0.0

    ra, dec = g("ra"), g("dec")
    if isinstance(ra, str):  # старые версии astroquery отдают sexagesimal
        c = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
        ra, dec = c.ra.deg, c.dec.deg
    return Target(id="", name=name, ra=float(ra), dec=float(dec),
                  pmra=g("pmra"), pmdec=g("pmdec"), parallax=g("plx_value"), rv=0.0,
                  ref_epoch=Time(2000.0, format="jyear", scale="tt"), source="SIMBAD")
