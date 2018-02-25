import datetime

a = datetime.datetime.strptime('02-26-02','%m-%d-%y')
print type(a)
print a.date()
# Raises a ValueError if format that doesn't match
a = datetime.datetime.strptime('02-26-02','%d/%m/%Y')
