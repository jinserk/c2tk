from pathlib import Path

import gpaw

from .. import settings


class GPAW(gpaw.GPAW):

    def __init__(self, *args, **kwargs):
        if 'directory' not in kwargs:
            kwargs['directory'] = settings.SCRATCH_PATH
        if 'txt' in kwargs:
            kwargs['txt'] = str(Path(kwargs['directory'], kwargs['txt']))
        super().__init__(*args, **kwargs)

    def write(self, filename, mode=''):
        filename = str(Path(self.directory, filename))
        super().write(filename, mode)

