import datetime
import json
import urllib
import urllib2

trends_url = 'http://api.twitter.com/1/trends/current.json'
last_week = datetime.date.today() - datetime.timedelta(days=7)
params = urllib.urlencode({
    'date': last_week.isoformat(),
    'exclude': 'hashtags'
})
url = '%s?%s' % (trends_url, params)
try:
    response = urllib2.urlopen(url, timeout=10)
except urllib2.HTTPError as e:
    print u'HTTPError getting Twitter data: %s' % e
except urllib2.URLError as e:
    print u'URLError getting Twitter data: %s' % e
else:
    content = response.read()
    data = json.loads(content)
    trends = data['trends']
    for trend_date, trend_list in trends.items():
        print u'Trends for %s' % trend_date
        for trend in trend_list:
            print trend['query']
