def fibonacci(n):
    if hasattr(fibonacci, '_%s' % n):
        return getattr(fibonacci, '_%s' % n)
    if n <= 2:
        value = 1
        setattr(fibonacci, '_%s' % n, value)
        return value
    else:
        value =  fibonacci(n - 1) + fibonacci(n - 2)
        setattr(fibonacci, '_%s' % n, value)
        return value

print fibonacci(10)
print fibonacci._10
print fibonacci(11)
print fibonacci._11
