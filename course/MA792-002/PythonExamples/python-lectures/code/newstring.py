class LanguageString(str):

    def __new__(cls, value, lang=u'en'):
        obj = str.__new__(cls, value)
        obj.lang = lang
        return obj

english_string = LanguageString('Hello')
spanish_string = LanguageString('Hola', lang='sp')
