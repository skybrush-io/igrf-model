from __future__ import annotations

from datetime import datetime
from importlib.resources import open_text
from pathlib import Path
from typing import ClassVar, Dict, IO, List, Optional, Tuple, Union

from ._igrf13syn import igrf13syn, InputType
from .vector import MagneticVector

IGRFModelCoefficients = Tuple[List[List[float]], List[List[float]]]

_DATA_PACKAGE = __name__.rpartition(".")[0] + ".data"


class IGRFModel:
    _cache: ClassVar[Dict[int, "IGRFModel"]] = {}

    @classmethod
    def from_file(cls, path: Union[str, Path]):
        with Path(path).open("r") as fp:
            return cls.from_io(fp)

    @classmethod
    def from_io(cls, fp: IO[str]):
        rows: List[List[float]] = []
        for line in fp:
            if line.startswith("g ") or line.startswith("h "):
                rows.append([float(x) for x in line.split()[3:]])

        transposed = list(map(list, zip(*rows)))

        coeffs = []
        for year_idx, row in enumerate(transposed):
            coeffs.extend(row[:120] if year_idx < 19 else row)

        coeffs.append(0)  # sentinel

        return cls(coeffs)

    @classmethod
    def get(cls, version: int = 13):
        instance = cls._cache.get(version)
        if instance is None:
            try:
                fp = open_text(_DATA_PACKAGE, f"igrf{version}.txt")
            except Exception as ex:
                raise ValueError("unsupported model version: {version}")
            with fp:
                instance = cls.from_io(fp)
            cls._cache[version] = instance
        return instance

    def __init__(
        self, table: List[float], *, min_year: float = 1900, max_year: float = 2030
    ):
        self._table = table
        self._min_year = float(min_year)
        # TODO(ntamas): calculate from self._table and self._min_year
        self._max_year = float(max_year)

    def evaluate(
        self,
        lat: float,
        lon: float,
        alt: float = 0.0,
        *,
        date: Optional[Union[float, datetime]] = None,
    ) -> MagneticVector:
        """Evaluates the model and calculates magnetic vector at the given latitude,
        longitude and altitude.

        Parameters:
            lat: the latitude (North is positive)
            lon: the longitude (East is positive)
            alt: the altitude above the WGS84 ellipsoid, in meters
            year: the (fractional) year to evaluate the model in; may also be a
                standard Python datetime object. `None` means to use the current
                date

        Returns:
            the magnetic vector
        """
        if date is None:
            date = datetime.now()
        return MagneticVector(
            *igrf13syn(
                date,
                InputType.GEODETIC,
                alt / 1000.0,
                lat,
                lon,
                gh=self._get_coeffs_for_year,
            )
        )

    def _get_coeffs_for_year(self, date: float) -> IGRFModelCoefficients:
        """Returns the coefficients of the IGRF magnetic model in the given year.

        Parameters:
            date: the year number; fractions are allowed. Valid range is between
                1900 and 2030, although accuracy is reduced beyond 2025.

        Returns:
            the G and H coefficients as two lists.
        """
        if date < self._min_year or date > self._max_year:
            raise ValueError(
                f"year out of range; allowed range is [{self._min_year}; {self._max_year}]"
            )

        date_threshold = self._max_year - 10

        # SH models before 1995.0 are only to degree 10
        nmx = 10 if date < 1995.0 else 13
        nc = nmx * (nmx + 2)

        if date >= date_threshold:
            # Extrapolation from date_threshold to the future
            t = date - date_threshold
            tc = 1.0
            ll = int((date_threshold - 1995) / 5)
            # 19 is the number of SH models that extend to degree 10
            ll = 120 * 19 + nc * ll
            # assertion: ll = 3255 if self._max_year == 2030, check it!
        else:
            # Interpolation between years
            t = (date - self._min_year) / 5
            ll = int(t)
            t = t - ll

            # SH models before 1995.0 are only to degree 10
            if date < 1995.0:
                ll = nc * ll
            else:
                ll = int((date - 1995) / 5)
                # 19 is the number of SH models that extend to degree 10
                ll = 120 * 19 + nc * ll
            tc = 1.0 - t

        table = self._table

        g, h = [], []
        temp = ll - 1
        for n in range(nmx + 1):
            g.append([])
            h.append([])
            if n == 0:
                g[0].append(None)
            for m in range(n + 1):
                if m != 0:
                    g[n].append(tc * table[temp] + t * table[temp + nc])
                    h[n].append(tc * table[temp + 1] + t * table[temp + nc + 1])
                    temp += 2
                else:
                    g[n].append(tc * table[temp] + t * table[temp + nc])
                    h[n].append(None)
                    temp += 1

        return g, h
