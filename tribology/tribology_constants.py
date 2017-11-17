# -*- coding: utf-8 -*-

"""

Relevant constants for tribology research

"""

from enum import Enum


class RadBall(Enum):
    """

    Constants for ball radii

    """
    TQInch = 25.4 * 3 / 8
    HInch = 25.4 / 4
    Inch = 25.4 / 2


class YoungsMod(Enum):
    """

    Young's modulus data in MPa at ambient temperature

    """
    STEEL = 210000
    GLASS = 70000
    SiN = 315


class PoissonRatio(Enum):
    """

    Poisson ratio data at ambient temperature

    """
    STEEL = 0.3
    GLASS = 0.22
    SiN = 0.26


class MatDens(Enum):
    """

    Density in kg / m$^3$ at ambient temperature

    """
    STEEL = 7800
    SiN = 3200


class LubeDens(Enum):
    """

    Density in g / ml

    """
    NA_LUBE_KR_015 = 0.884
    NYNAS_HP_4 = 0.852
    NYNAS_HP_12 = 0.867
    NYNAS_T_3 = 0.868
    NYNAS_T_9 = 0.888
    NYNAS_T_22 = 0.902
    SIGMA_ALDRICH_MINERAL_OIL_HEAVY = 0.862


class LubeVisc(Enum):
    """

    Lubricant viscosity in cSt at 40 and 100 degree C

    """
    TEMPS = (40, 100)
    NA_LUBE_KR_015 = (114.0, 13.5)
    NYNAS_HP_4 = (20, 4.2)
    NYNAS_HP_12 = (110, 12)
    NYNAS_T_3 = (3.7, 1.3)
    NYNAS_T_9 = (9.0, 2.2)
    NYNAS_T_22 = (22.5, 3.6)
    SIGMA_ALDRICH_MINERAL_OIL_HEAVY = (67.0, 18.9)


class PressVisc(Enum):
    """

    Pressure-viscosity coefficients in 1 / Pa at ambient temperature

    """
    ESTER_OIL_GENERIC = 15 * 10 ** (-9)
    MINERAL_OIL_GENERIC = 30 * 10 ** (-9)
