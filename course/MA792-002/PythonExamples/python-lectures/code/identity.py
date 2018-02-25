def identity(n):
    return [[1 if i == j else 0 for i in xrange(n)] for j in xrange(n)]

print identity(3)

print [(i, j) for i in xrange(3) for j in xrange(3)]
