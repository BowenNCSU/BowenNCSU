def simple(x, y=[]):
    y.append(x)
    return y

print simple.func_defaults
print simple(1)
print simple(2)
print simple.func_defaults
