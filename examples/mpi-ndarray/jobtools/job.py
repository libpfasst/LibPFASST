
from base import Container, NoInputError


class Job(Container):
    """Job container."""

    def update_params(self, **kwargs):
        if self.param_file is not None:
            with open(self.param_file, 'r') as f:
                template = f.read()
            self._inputs = template.format(**kwargs)
        else:
            raise NoInputError('update_params called but param_file not set.')


    def write_params(self, fname):
        if self._inputs is not None:
            with open(fname, 'w') as f:
                f.write(self._inputs)

    @property
    def has_params(self):
        return self._inputs is not None
