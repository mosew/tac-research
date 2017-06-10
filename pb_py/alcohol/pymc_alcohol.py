from pymc import *

q1mean = Uniform('q1mean', lower = .00001, upper = .001)
q1precision = Uniform('q1precision', lower = .00001, upper = .001)
