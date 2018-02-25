class Student(object):

    def __init__(self, name, grade):
        self.name, self.grade = name, grade

    def __lt__(self, other):
        return self.name < other.name

    def __cmp__(self, other):
        return self.grade - other.grade


sally = Student('Sally', 90)
susan = Student('Susan', 90)
print sally < susan
print sally > susan
print sally >= susan
