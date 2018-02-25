class Person(object):

    def __init__(self, name, age):
        self.name = name
        self.age = age

    def is_old(self):
        return self.age > 40

        
person = Person('G. H. Hardy', 70)
print person.is_old()
