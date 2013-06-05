
class Container(object):
    def __init__(self, **kwargs):
        self.attrs = {}
        self.attrs.update(kwargs)
    def __getattr__(self, name):
        return self.attrs.get(name, None)

class NoInputError(Exception):
    pass
