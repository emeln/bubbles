from .bubble_functions import BubbleDist
from .sigma_interpolation import SigmaInterpolator
from .cosmology import Cosmology

bubblyTest = BubbleDist()
sigmaInt = SigmaInterpolator()
m1 = 10e-24
z1 = 1
sigmaInt(m1, z1)