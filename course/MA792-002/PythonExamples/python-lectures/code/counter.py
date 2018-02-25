class Counter(object):
    count = 0                     
    def __init__(self):
        self.__class__.count += 1

first = Counter()
print first.count
second = Counter()
print second.count
