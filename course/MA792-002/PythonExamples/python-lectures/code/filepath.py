import os

print __file__
print os.getcwd()
print os.path.abspath(__file__)
print os.path.abspath(os.path.join(os.getcwd(), '..', 'test'))
print os.path.splitext('test.txt')
