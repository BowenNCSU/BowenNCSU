def greater_than_bound(bound):
    def greater(x):
        return x > bound
    return greater

greater_than_ten = greater_than_bound(10)
print greater_than_ten(20)

