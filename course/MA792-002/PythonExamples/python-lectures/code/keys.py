from operator import itemgetter, attrgetter, methodcaller

test = [(1, 2), (5, 3), (7, 1), (3, 8)]
print sorted(test, key=itemgetter(1))

test = ['C', 'a', 'd', 'B', 'E']
print sorted(test, key=methodcaller('lower'))

test = [1 +  1j, 2 - 1j, -1 + 2j]
print sorted(test, key=attrgetter('real'))
