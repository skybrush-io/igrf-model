from datetime import datetime
from typing import Optional, Union

__all__ = ("FractionalYearLike",)

FractionalYearLike = Optional[Union[datetime, float]]
"""Type alias for objects that can be converted into a fractional year."""
