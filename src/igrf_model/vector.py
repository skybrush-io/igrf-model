from dataclasses import dataclass
from math import atan2, degrees, hypot, sqrt

__all__ = ("MagneticVector",)


@dataclass
class MagneticVector:
    """Simple dataclass representing a magnetic vector."""

    north: float
    """The North component of the vector."""

    east: float
    """The East component of the vector."""

    down: float
    """The vertical component of the vector (positive points down)."""

    @property
    def declination(self) -> float:
        """The declination of the magnetic vector, in degrees (positive points East)."""
        return degrees(atan2(self.east, self.north))

    @property
    def inclination(self) -> float:
        """The inclination of the magnetic vector, in degrees (positive points down)."""
        return degrees(atan2(self.down, self.horizontal_intensity))

    @property
    def horizontal_intensity(self) -> float:
        """The horizontal intensity of the magnetic vector."""
        return hypot(self.north, self.east)

    @property
    def total_intensity(self) -> float:
        """The total intensity of the magnetic vector."""
        return sqrt(
            self.north * self.north + self.east * self.east + self.down * self.down
        )
