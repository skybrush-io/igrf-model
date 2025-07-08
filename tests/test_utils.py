from datetime import datetime
from pytest import approx

from igrf_model.utils import (
    datetime_to_fractional_year,
    parse_datetime_or_fractional_year,
)


def test_datetime_to_fractional_year():
    d2y = datetime_to_fractional_year
    assert d2y(datetime(2010, 1, 1)) == 2010
    assert d2y(datetime(2007, 12, 31)) == 2007 + 364 / 365
    assert d2y(datetime(2017, 4, 1)) == 2017 + 90 / 365


def test_parse_datetime_or_fractional_year():
    now = datetime.now()
    parse = parse_datetime_or_fractional_year
    assert parse(None) == approx(parse(now))
    assert parse(1974.5) == 1974.5
    assert parse(datetime(2017, 4, 1)) == 2017 + 90 / 365
