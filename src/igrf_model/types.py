from datetime import datetime
from typing import List, Optional, Tuple, Union

__all__ = ("FractionalYearLike",)

FractionalYearLike = Optional[Union[datetime, float]]
"""Type alias for objects that can be converted into a fractional year."""

IGRFModelCoefficients = Tuple[List[List[float]], List[List[float]]]
"""Type alias for the IGRF model coefficients (g and h)"""
