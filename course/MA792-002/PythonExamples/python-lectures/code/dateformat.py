import datetime

today = datetime.date.today()
now = datetime.datetime.now()
print today.strftime('%a %B %d, %Y')
print now.strftime('%m-%d-%y %I:%M %p')
