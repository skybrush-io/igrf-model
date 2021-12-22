from .model import DateBoundIGRFModel, IGRFModel
from .types import FractionalYearLike
from .vector import MagneticVector
from .version import __version__, __version_info__

__all__ = (
    "DateBoundIGRFModel",
    "IGRFModel",
    "MagneticVector",
    "FractionalYearLike",
    "__version__",
    "__version_info__",
)
