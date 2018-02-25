from itertools import chain, takewhile, dropwhile
from itertools import combinations, permutations

a = [1, 2, 3]
b = 'abc'
print [(x, type(x)) for x in chain(a, b)]
print [u''.join(x) for x in combinations(b, 2)]
print [u''.join(x) for x in permutations(b, 2)]
print list(takewhile(lambda x: x % 2 == 1, a))
print list(dropwhile(lambda x: x in 'aeiou', b))
