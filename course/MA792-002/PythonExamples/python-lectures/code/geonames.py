"""
Parsing xml data from GeoNames
http://www.geonames.org/export/web-services.html
"""

import os
import urllib2

from lxml import etree


def save_country_info(filename):
    url = 'http://api.geonames.org/countryInfo?username=demo'
    try:
        response = urllib2.urlopen(url, timeout=10)
    except urllib2.HTTPError as e:
        print u'HTTPError getting GeoNames data: %s' % e
    except urllib2.URLError as e:
        print u'URLError getting GeoNames data: %s' % e
    else:
        with open(filename, 'w') as f:
            f.write(response.read())


def parse_country_info(filename):
    with open(filename) as f:
        tree = etree.parse(f)
    return tree


if __name__ == '__main__':
    cache_file = 'geonames.xml'
    if not os.path.exists(cache_file):
        save_country_info(cache_file)
    tree = parse_country_info(cache_file)
    root = tree.getroot()
    countries = root.xpath('//country[population>1000000]')
    for country in countries:
        name = country.findtext('countryName')
        population = country.findtext('population')
        print u'%s - Population %d' % (name, int(population))
