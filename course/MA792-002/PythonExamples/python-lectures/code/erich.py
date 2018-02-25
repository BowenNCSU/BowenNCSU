class Example(object):

    def test(self):
        print 'blip'

def new_test():
    print 'new'

example = Example()
example.test()

example.test = new_test
example.test()
example = Example()
example.test()


