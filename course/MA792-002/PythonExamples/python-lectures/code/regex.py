import random
import re
from string import ascii_lowercase, digits

simple_phone = re.compile(r'[1-9]\d{2}(-|\.|\s)\d{4}')
garbage = ''.join([
    random.choice(ascii_lowercase + digits + '-. ')
    for x in range(20000)
])
match = simple_phone.search(garbage)
if match:
    print garbage
    print 'Found a match!'
    print garbage[match.start():match.end()]
else:
    print 'No matches...'
