# -*- coding: utf-8 -*-

"""

Relevant constants for tribology research. Use, for example, as follows:

.. code-block:: python

   # get young's modulus of steel in units of MPa
   import tribology as tr
   youngs_modulus_steel = tr.YoungsMod.STEEL.value

"""

from enum import Enum


class RadBall(Enum):
    """

    Scalar constants for ball radii in units of mm. Available values are:

    :TQInch: 3/4 inch
    :HInch: 1/2 inch
    :Inch: 1 inch

    """
    TQInch = 25.4 * 3 / 8
    HInch = 25.4 / 4
    Inch = 25.4 / 2


class YoungsMod(Enum):
    """

    Scalar constants for Young's modulus in units of MPa. Available values are:

    :STEEL:
    :GLASS:
    :SiN:

    """
    STEEL = 210000
    GLASS = 70000
    SiN = 315000


class PoissonRatio(Enum):
    """

    Scalar constants for unitless Poisson ratios. Available values are:

    :STEEL:
    :GLASS:
    :SiN:

    """
    STEEL = 0.3
    GLASS = 0.22
    SiN = 0.26


class MatDens(Enum):
    """

    Scalar constants for material density in :math:`\\text{kg m}^{-3}`.
    Available values are:

    :STEEL:
    :GLASS:

    """
    STEEL = 7800
    SiN = 3200


class LubeDens(Enum):
    """

    Scalar constants for lubricant density in :math:`\\text{g ml}^{-1}`.
    Available values are:

    :NA_LUBE_KR_015:
    :NYNAS_HP_4:
    :NYNAS_HP_12:
    :NYNAS_T_3:
    :NYNAS_T_9:
    :NYNAS_T_22:
    :SIGMA_ALDRICH_MINERAL_OIL_HEAVY:

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

    Tuple containing kinematic viscosities (in unit of cSt) of lubricants
    at 40 and 100 :math:`^{\\circ}\\text{C}`. The first value in the tuple
    corresponds to the viscosity at 40 :math:`^{\\circ}\\text{C}`, the second
    to 100 :math:`^{\\circ}\\text{C}`. The **TEMPS** value contains the
    reference temperatures (40, 100). Available values are:

    :TEMPS: Tuple containing temperature values of 40 and 100

    :NA_LUBE_KR_015:
    :NYNAS_HP_4:
    :NYNAS_HP_12:
    :NYNAS_T_3:
    :NYNAS_T_9:
    :NYNAS_T_22:
    :SIGMA_ALDRICH_MINERAL_OIL_HEAVY:

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

    Pressure-viscosity coefficients in :math:`\\text{Pa}^{-1}` at room
    temperature. Available values are:

    :ESTER_OIL_GENERIC:
    :MINERAL_OIL_GENERIC:

    """
    ESTER_OIL_GENERIC = 15 * 10 ** (-9)
    MINERAL_OIL_GENERIC = 30 * 10 ** (-9)
