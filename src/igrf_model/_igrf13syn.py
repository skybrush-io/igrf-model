from __future__ import annotations

from datetime import date, datetime
from enum import IntEnum
from math import atan2, cos, radians, sin, sqrt
from typing import Callable, Tuple, Union, TYPE_CHECKING

if TYPE_CHECKING:
    from .model import IGRFModelCoefficients  # pragma: no cover


class InputType(IntEnum):
    """Enum representation of the possible input types of the `igrf13syn`
    function, converted from Fortran.
    """

    GEODETIC = 1
    """Geodetic coordinates (WGS84 spheroid)."""

    GEOCENTRIC = 2
    """Geocentric coordinates."""


def datetime_to_fractional_year(input: datetime) -> float:
    """Converts a Python datetime object to a fractional year."""
    start = date(input.year, 1, 1).toordinal()  # type: ignore
    year_length = date(input.year + 1, 1, 1).toordinal() - start  # type: ignore
    return input.year + (input.toordinal() - start) / year_length


def sin_and_cos(radians: float) -> Tuple[float, float]:
    """Helper function to compute the sine and cosine of an angle simultaneously."""
    return sin(radians), cos(radians)


def igrf13syn(
    date: Union[datetime, float],
    itype: InputType,
    alt: float,
    lat: float,
    elong: float,
    *,
    gh: Callable[[float], "IGRFModelCoefficients"]
) -> Tuple[float, float, float]:
    """This is a synthesis routine for the 13th generation IGRF as agreed
    in December 2019 by IAGA Working Group V-MOD. It is valid 1900.0 to
    2025.0 inclusive. Values for dates from 1945.0 to 2015.0 inclusive are
    definitive, otherwise they are non-definitive.

    Parameters:
        date: year A.D. Must be greater than or equal to 1900.0 and
            less than or equal to 2030.0. May be fractional. You can also
            provide a standard Python datetime object here; it will be converted
            appropriately.
        itype: whether the input coordinates are geodetic or geocentric. Use
            the constants from the InputType enum for clarity.
        alt: height in km above sea level if itype = 1; distance from centre of
            Earth in km if itype = 2 (>3485 km)
        lat: latitude (-90 to 90)
        elong: east-longitude (0-360)
        gh: a function that can be called with a fractional year and that
            returns the IGRF model coefficients

    Returns:
        a tuple consisting of the North, East and vertical components (+ve down)
        of the magnetic vector (in nT)

    This function was converted from the Fortran IGRF13SYN function from
    https://github.com/space-physics/igrf/blob/main/src/igrf/fortran/igrf13.f .
    Its interface was modified slightly to be more Pythonic. The comments
    below are from the original Fortran implementation:

    Adapted from 8th generation version to include new maximum degree for
    main-field models for 2000.0 and onwards and use WGS84 spheroid instead
    of International Astronomical Union 1966 spheroid as recommended by IAGA
    in July 2003. Reference radius remains as 6371.2 km - it is NOT the mean
    radius (= 6371.0 km) but 6371.2 km is what is used in determining the
    coefficients. Adaptation by Susan Macmillan, August 2003 (for
    9th generation), December 2004, December 2009, December 2014;
    by William Brown, December 2019, February 2020.

    Coefficients at 1995.0 incorrectly rounded (rounded up instead of
    to even) included as these are the coefficients published in Excel
    spreadsheet July 2005.
    """
    if not isinstance(date, (int, float)):
        date = datetime_to_fractional_year(date)
    else:
        date = float(date)

    g, h = gh(date)

    p, q, cl, sl = [0.0] * 105, [0.0] * 105, [0.0] * 13, [0.0] * 13
    x, y, z = 0.0, 0.0, 0.0

    nmx = len(g) - 1
    kmx = (nmx + 1) * (nmx + 2) // 2

    # Conversion from our parameters to the ones used by the original Fortran
    # routine
    colat = 90 - lat
    r = alt

    # Fortran uses colat*0.017453292 but this is more readable
    st, ct = sin_and_cos(radians(colat))
    sl[0], cl[0] = sin_and_cos(radians(elong))
    l, m, n = 1, 1, 0

    if itype is InputType.GEODETIC:
        # conversion from geodetic to geocentric coordinates
        # (using the WGS84 spheroid)
        a2 = 40680631.6
        b2 = 40408296.0
        one = a2 * st * st
        two = b2 * ct * ct
        three = one + two
        rho = sqrt(three)
        r = sqrt(alt * (alt + 2.0 * rho) + (a2 * one + b2 * two) / three)
        cd = (alt + rho) / r
        sd = (a2 - b2) / rho * ct * st / r
        one = ct
        ct = ct * cd - st * sd
        st = st * cd + one * sd
    else:
        cd, sd = 1.0, 0.0

    ratio = 6371.2 / r
    rr = ratio * ratio

    # computation of Schmidt quasi-normal coefficients p and x(=q)
    p[0] = 1.0
    p[2] = st
    q[0] = 0.0
    q[2] = ct

    # fn, gn lifted here from the loop because they never change
    fn, gn = n, n - 1

    for k in range(2, kmx + 1):
        if n < m:
            m = 0
            n = n + 1
            rr = rr * ratio
            fn = n
            gn = n - 1

        fm = m
        if m != n:
            gmm = m * m
            one = sqrt(fn * fn - gmm)
            two = sqrt(gn * gn - gmm) / one
            three = (fn + gn) / one
            i = k - n
            j = i - n + 1
            p[k - 1] = three * ct * p[i - 1] - two * p[j - 1]
            q[k - 1] = three * (ct * q[i - 1] - st * p[i - 1]) - two * q[j - 1]
        elif k != 3:
            one = sqrt(1.0 - 0.5 / fm)
            j = k - n - 1
            p[k - 1] = one * st * p[j - 1]
            q[k - 1] = one * (st * q[j - 1] + ct * p[j - 1])
            cl[m - 1] = cl[m - 2] * cl[0] - sl[m - 2] * sl[0]
            sl[m - 1] = sl[m - 2] * cl[0] + cl[m - 2] * sl[0]

        # synthesis of x, y and z in geocentric coordinates
        one = g[n][m] * rr
        if m == 0:
            x = x + one * q[k - 1]
            z = z - (fn + 1.0) * one * p[k - 1]
            l = l + 1
        else:
            two = h[n][m] * rr
            three = one * cl[m - 1] + two * sl[m - 1]
            x = x + three * q[k - 1]
            z = z - (fn + 1.0) * three * p[k - 1]
            if st == 0:
                y = y + (one * sl[m - 1] - two * cl[m - 1]) * q[k - 1] * ct
            else:
                y = y + (one * sl[m - 1] - two * cl[m - 1]) * fm * p[k - 1] / st
            l = l + 2
        m = m + 1

    # conversion to coordinate system specified by itype
    one = x
    x = x * cd + z * sd
    z = z * cd - one * sd
    return x, y, z
