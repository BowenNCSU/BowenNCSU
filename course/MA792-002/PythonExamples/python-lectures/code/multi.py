class TypeA(object):

    def name(self):
        print u"Type A"


class TypeB(object):

    def name(self):
        print u"Type B"


class TypeC(TypeA, TypeB):
    pass


class TypeD(TypeB, TypeA):
    pass


c = TypeC()
c.name()
print c.__class__.__mro__

d = TypeD()
d.name()
print d.__class__.__mro__

