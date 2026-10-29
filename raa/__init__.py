"""Revised affine arithmetic (RAA) and the AF2 baseline.

RAA represents every subexpression as a linear function of the input variables
plus a bounded residual interval, and re-centers the residuals at every
multiplication. See the README for the reference to the paper.
"""
from .box import set_box, get_box
from .raa import RAAForm, raa_variables, pow_raa, mul_raa
from .af2 import AF2, af2_variables
from .unary import sin, cos, exp, log, sqrt, sinh, power

__all__ = ["set_box", "get_box", "RAAForm", "raa_variables", "pow_raa", "mul_raa",
           "AF2", "af2_variables", "sin", "cos", "exp", "log", "sqrt", "sinh", "power"]
