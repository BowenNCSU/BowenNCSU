from scipy import array
from scipy.stats import linregress

from geonames import parse_country_info


def get_country_data(country):
    size = country.findtext('areaInSqKm')
    population = country.findtext('population')
    if size and population:
        return size, population
    else:
        return None, None

def parse_pop_size_data():
    tree = parse_country_info('geonames.xml')
    root = tree.getroot()
    data = map(get_country_data, root.findall('country'))
    data = array(filter(lambda x: x[0] and x[1], data), dtype=float)
    data.sort()
    size = data[:,0]
    population = data[:,1]
    return (size, population)


if __name__ == '__main__':
    (a_s, b_s, r, tt, stderr) = linregress(*parse_pop_size_data())
    print 'Regression: a=%.2f b=%.2f, r=%.2f, std error= %.3f' % (a_s, b_s, r, stderr)
