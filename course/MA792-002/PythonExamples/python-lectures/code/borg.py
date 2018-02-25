class Borg(object):
    _state = {}
    
    def __new__(cls, *p, **k):
        self = object.__new__(cls, *p, **k)
        # override instance namespace with shared state
        self.__dict__ = cls._state
        return self


config1 = Borg()
config1.debug = True

config2 = Borg()
print config2.debug

print config1 is config2
