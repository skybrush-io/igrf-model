from math import degrees, hypot

from ._igrf13syn import igrf13syn, InputType


def igrf_variation(lat, lon, alt=0.0, year=2005):
    """
    Annual variation
    D is declination (+ve east)
    I is inclination (+ve down)
    H is horizontal intensity
    x is north component
    y is east component
    Z is vertical component (+ve down)
    F is total intensity
    """
    x1, y1, z1, f1 = igrf13syn(year - 1, InputType.GEODETIC, alt, lat, lon)
    x2, y2, z2, f2 = igrf13syn(year + 1, InputType.GEODETIC, alt, lat, lon)
    x, y, z, f = (x1 + x2) / 2, (y1 + y2) / 2, (z1 + z2) / 2, (f1 + f2) / 2
    dx, dy, dz, df = (x2 - x1) / 2, (y2 - y1) / 2, (z2 - z1) / 2, (f2 - f1) / 2
    h = hypot(x, y)

    dd = degrees(x * dy - y * dx) / (h * h)
    dh = (x * dx + y * dy) / h
    ds = degrees(h * dz - z * dh) / (f * f)
    df = (h * dh + z * dz) / f
    return dd, ds, dh, dx, dy, dz, df
