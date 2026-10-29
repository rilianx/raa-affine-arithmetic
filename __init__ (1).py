"""Zero-centered affine arithmetic (ZCAA) and the AF2 and first-order Taylor model baselines.

ZCAA represents every subexpression as a linear function of the input variables
plus a bounded residual interval, and re-centers the residuals at every
multiplication. See the README for the reference to the paper.
"""
from .box import set_box, get_box
from .zcaa import ZCAAForm, zcaa_variables, pow_zcaa, mul_zcaa
from .af2 import AF2, af2_variables
from .tm import TM1, tm_variables
from .unary import sin, cos, exp, log, sqrt, sinh, power

__all__ = ["set_box", "get_box", "ZCAAForm", "zcaa_variables", "pow_zcaa", "mul_zcaa",
           "AF2", "af2_variables", "TM1", "tm_variables", "sin", "cos", "exp", "log", "sqrt", "sinh", "power"]
