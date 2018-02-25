import os

root = os.getcwd()
for item in os.listdir(root):
    print item

def visit(arg, dirname, names):
    print arg, dirname, names

os.path.walk(root, visit, None)
