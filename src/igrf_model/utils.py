from datetime import date, datetime

from .types import FractionalYearLike

__all__ = ("datetime_to_fractional_year", "parse_datetime_or_fractional_year")


def datetime_to_fractional_year(input: datetime) -> float:
    """Converts a Python datetime object to a fractional year."""
    start = date(input.year, 1, 1).toordinal()  # type: ignore
    year_length = date(input.year + 1, 1, 1).toordinal() - start  # type: ignore
    return input.year + (input.toordinal() - start) / year_length


def parse_datetime_or_fractional_year(input: FractionalYearLike) -> float:
    """Converts a Python datetime object to a fractional year; also accepts
    `None` to represent the current date and any single floating-point number
    that is interpreted as a fractional year as is.
    """
    if input is None:
        input = datetime.now()
    if isinstance(input, datetime):
        input = datetime_to_fractional_year(input)
    return input
