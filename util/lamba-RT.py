import math
import numpy
from scipy.special import lambertw

# parameters for R
increment = 0.005
limit = 1.5
base = 1 + increment
current = base

# increment the current value of R until the limit is reached
while current < limit:
    argument = - current / math.exp(current)
    lambert = numpy.real(lambertw(argument))

    # the exacti value for lambda * T
    value = current + lambert

    ## estimations for lambda * T

    # this estimation was empirical at first but it is also an agressive simplification of the estimation below
    estimation  = (current-1)*2

    # this estimation comes from https://link.springer.com/article/10.1007/s10444-017-9530-3 , equation 3
    estimation2 = current - 1 + math.sqrt(2)*math.sqrt(1 + math.exp(1)*argument)

    ratio  = estimation / value
    ratio2 = estimation2 / value

    print( current, estimation, estimation2, lambert, value, ratio, ratio2 )

    current = current + increment

