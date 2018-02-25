size = 10
print range(size)
L = range(size)
for i in range(0,size):
   print L[i]

from array import array
arr = array('i')
for i in range(size):
    arr.append(L[i])
print arr

text = '1,0,0;0,1,0;0,0,1'
rows = text.split(';')
matrix = []
for row in rows:
    matrix.append([float(x) for x in row.split(',')])
print matrix

matrix = map(lambda x: map(float, x.split(',')), text.split(';'))
print matrix
