
from fabric.utils import _AttributeDict

class Container(_AttributeDict):
    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError:
            return None

class NoInputError(Exception):
    pass
