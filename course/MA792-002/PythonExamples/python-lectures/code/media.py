class Media(object):
    
    def __init__(self, css=None, js=None):
        self.css = set(css or [])
        self.js = set(js or [])

    def __add__(self, other):
        css = self.css | other.css
        js = self.js | other.js
        return Media(css=css, js=js)

base = Media(css=['base.css'])
forms = Media(css=['base.css', 'forms.css'], js=['forms.js'])

new = base + forms
print new.css
print new.js
