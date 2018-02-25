from person import Person

class Student(Person):
    
    def __init__(self, name, age, gpa):
        self.gpa = gpa
        super(Student, self).__init__(name, age)

    def is_honor_student(self):
        return self.gpa > 3.0


student = Student('G. H. Hardy', 70, 4.0)
print student.is_old()
print student.is_honor_student()
