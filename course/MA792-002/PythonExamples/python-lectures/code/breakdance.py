for name in ['Newton', 'x', 'Euler']:
    if 'x' in name:
        continue
    else:
        print name[0]

for name in []:
    print "Names were given."
else:
    print "No names were given."

for name in ['Hilbert', 'Mark', 'Hardy']:
    if name == 'Mark':
        print "Mark is the greatest!"
        break
    else:
        print name[0]
