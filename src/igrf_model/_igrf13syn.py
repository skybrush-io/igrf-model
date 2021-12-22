from __future__ import annotations

from enum import IntEnum
from math import cos, radians, sin, sqrt
from typing import Callable, Tuple

from .types import FractionalYearLike, IGRFModelCoefficients
from .utils import parse_datetime_or_fractional_year

__all__ = ("igrf13syn", "igrf13syn_with_precalculated_coeffs", "InputType")


class InputType(IntEnum):
    """Enum representation of the possible input types of the `igrf13syn`
    function, converted from Fortran.
    """

    GEODETIC = 1
    """Geodetic coordinates (WGS84 spheroid)."""

    GEOCENTRIC = 2
    """Geocentric coordinates."""


def sin_and_cos(radians: float) -> Tuple[float, float]:
    """Helper function to compute the sine and cosine of an angle simultaneously."""
    return sin(radians), cos(radians)


def igrf13syn(
    date: FractionalYearLike,
    itype: InputType,
    alt: float,
    lat: float,
    elong: float,
    gh: Callable[[float], IGRFModelCoefficients],
) -> Tuple[float, float, float]:
    """This is a synthesis routine for the 13th generation IGRF as agreed
    in December 2019 by IAGA Working Group V-MOD. It is valid 1900.0 to
    2025.0 inclusive. Values for dates from 1945.0 to 2015.0 inclusive are
    definitive, otherwise they are non-definitive.

    Parameters:
        date: year A.D. Must be greater than or equal to 1900.0 and
            less than or equal to 2030.0. May be fractional. You can also
            provide a standard Python datetime object here; it will be converted
            appropriately. `None` means to use the current date and time.
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
    Its interface was modified slightly to be more Pythonic and to cater for the
    common use-case when multiple points are evaluated at the same date (that is
    why the function can take the IGRF model coefficients at a given date
    directly). The comments below are from the original Fortran implementation:

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
    coeffs = gh(parse_datetime_or_fractional_year(date))
    return igrf13syn_with_precalculated_coeffs(coeffs, itype, alt, lat, elong)


def igrf13syn_with_precalculated_coeffs(
    coeffs: IGRFModelCoefficients,
    itype: InputType,
    alt: float,
    lat: float,
    elong: float,
) -> Tuple[float, float, float]:
    """This is a synthesis routine for the 13th generation IGRF as agreed
    in December 2019 by IAGA Working Group V-MOD. It is valid 1900.0 to
    2025.0 inclusive. Values for dates from 1945.0 to 2015.0 inclusive are
    definitive, otherwise they are non-definitive.

    Parameters:
        coeffs: the IGRF model coefficients for the date that the caller is
            interested in. Use `igrf13syn()` instead if you have a date.
        itype: whether the input coordinates are geodetic or geocentric. Use
            the constants from the InputType enum for clarity.
        alt: height in km above sea level if itype = 1; distance from centre of
            Earth in km if itype = 2 (>3485 km)
        lat: latitude (-90 to 90)
        elong: east-longitude (0-360)

    Returns:
        a tuple consisting of the North, East and vertical components (+ve down)
        of the magnetic vector (in nT)

    This function was converted from the Fortran IGRF13SYN function from
    https://github.com/space-physics/igrf/blob/main/src/igrf/fortran/igrf13.f .
    Its interface was modified slightly to be more Pythonic and to cater for the
    common use-case when multiple points are evaluated at the same date (that is
    why the function can take the IGRF model coefficients at a given date
    directly). The comments below are from the original Fortran implementation:

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

    p, q, cl, sl = [0.0] * 105, [0.0] * 105, [0.0] * 13, [0.0] * 13
    x, y, z = 0.0, 0.0, 0.0

    g, h = coeffs

    nmx = len(g) - 1
    kmx = (nmx + 1) * (nmx + 2) // 2

    # Conversion from our parameters to the ones used by the original Fortran
    # routine
    colat = 90 - lat
    r = alt

    # Fortran uses colat*0.017453292 but this is more readable
    st, ct = sin_and_cos(radians(colat))
    sl[0], cl[0] = sin_and_cos(radians(elong))

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

    m, n = 1, 0

    for k in range(1, kmx):
        if n < m:
            m = 0
            n = n + 1
            rr = rr * ratio

        if m != n:
            gmm = m * m
            one = sqrt(n * n - gmm)
            two = sqrt((n - 1) * (n - 1) - gmm) / one
            three = (2 * n - 1) / one
            i = k - n
            j = i - n + 1
            p[k] = three * ct * p[i] - two * p[j]
            q[k] = three * (ct * q[i] - st * p[i]) - two * q[j]
        elif k != 2:
            one = sqrt(1.0 - 0.5 / m)
            j = k - n - 1
            p[k] = one * st * p[j]
            q[k] = one * (st * q[j] + ct * p[j])
            cl[m - 1] = cl[m - 2] * cl[0] - sl[m - 2] * sl[0]
            sl[m - 1] = sl[m - 2] * cl[0] + cl[m - 2] * sl[0]

        # synthesis of x, y and z in geocentric coordinates
        one = g[n][m] * rr
        if m == 0:
            x = x + one * q[k]
            z = z - (n + 1) * one * p[k]
        else:
            two = h[n][m] * rr
            three = one * cl[m - 1] + two * sl[m - 1]
            x = x + three * q[k]
            z = z - (n + 1) * three * p[k]
            if st == 0:
                y = y + (one * sl[m - 1] - two * cl[m - 1]) * q[k] * ct
            else:
                y = y + (one * sl[m - 1] - two * cl[m - 1]) * m * p[k] / st

        m = m + 1

    # conversion to coordinate system specified by itype
    one = x
    x = x * cd + z * sd
    z = z * cd - one * sd
    return x, y, z
