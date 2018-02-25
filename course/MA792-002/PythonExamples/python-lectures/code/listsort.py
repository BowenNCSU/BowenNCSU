test = ['C', 'a', 'd', 'B', 'E']
test.sort()
print test # Default alpha sort

test.sort(key=str.lower)
print test # Case insensitive alpha sort

test = [(1, 2), (5, 3), (7, 1), (3, 8)]
test.sort(key=lambda x: x[1])
print test # Sort by second entry
