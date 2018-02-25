from matplotlib import pylab
from scipy import polyval
from scipy.stats import linregress

from regress import parse_pop_size_data


size, population = parse_pop_size_data()
(a_s, b_s, r, tt, stderr) = linregress(size, population)
line = polyval([a_s, b_s], size)

pylab.title('Linear Regression Example')
pylab.plot(size, population, 'r.', size, line, 'k')
pylab.xlabel('Area in sqKm')
pylab.ylabel('Population')
pylab.legend(['data', 'regression'])
pylab.show()
